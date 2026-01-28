#!/bin/python

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

import torch
import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.ops.indexing import Vindex
from pyro.infer import config_enumerate, infer_discrete
from pyro.infer.autoguide import AutoDiagonalNormal, AutoDelta
from pyro.infer import JitTraceEnum_ELBO, TraceEnum_ELBO
from pyro.infer.svi import SVI
from pyro.infer import MCMC, NUTS, HMC


class AneuploidyData:
    """Data container for aneuploidy detection model"""
    def __init__(self, df: pd.DataFrame, device: str = 'cpu', dtype: torch.dtype = torch.float32,
                 subsample_bins: int = None, subsample_samples: int = None, seed: int = 42,
                 clamp_threshold: float = 5.0):
        """
        Args:
            df: DataFrame with bins as rows and samples as columns
            device: torch device
            dtype: torch data type
            subsample_bins: If specified, randomly subsample this many bins
            subsample_samples: If specified, randomly subsample this many samples
            seed: Random seed for subsampling
            clamp_threshold: Maximum value for depth; values above this are clamped
        """
        # Store original dataframe
        self.original_df = df
        
        # Get sample columns (exclude metadata)
        metadata_cols = ['Chr', 'Start', 'End', 'source_file']
        sample_cols = [col for col in df.columns if col not in metadata_cols]
        
        # Subsample if requested
        if subsample_bins is not None or subsample_samples is not None:
            np.random.seed(seed)
            
            # Subsample bins (rows)
            if subsample_bins is not None and subsample_bins < len(df):
                print(f"Subsampling {subsample_bins} bins from {len(df)} total bins")
                bin_indices = np.random.choice(len(df), size=subsample_bins, replace=False)
                bin_indices = np.sort(bin_indices)  # Keep sorted for interpretability
                df = df.iloc[bin_indices].copy()
            
            # Subsample samples (columns)
            if subsample_samples is not None and subsample_samples < len(sample_cols):
                print(f"Subsampling {subsample_samples} samples from {len(sample_cols)} total samples")
                sample_indices = np.random.choice(len(sample_cols), size=subsample_samples, replace=False)
                selected_samples = [sample_cols[i] for i in sample_indices]
                sample_cols = selected_samples
        
        # Sort chromosomes in proper order (chr1-chr22, chrX, chrY)
        chr_order = {f'chr{i}': i for i in range(1, 23)}
        chr_order['chrX'] = 23
        chr_order['chrY'] = 24
        df['_chr_order'] = df['Chr'].map(chr_order)
        df = df.sort_values(['_chr_order', 'Start']).drop('_chr_order', axis=1)
        
        # Extract metadata
        self.chr = df['Chr'].values
        self.start = df['Start'].values
        self.end = df['End'].values
        self.sample_ids = sample_cols
        
        # Extract normalized read depth matrix (bins x samples)
        depth_matrix = df[sample_cols].values
        
        # Clamp values above threshold
        if clamp_threshold is not None:
            n_clamped = np.sum(depth_matrix > clamp_threshold)
            if n_clamped > 0:
                print(f"Clamping {n_clamped} values above threshold {clamp_threshold}")
                depth_matrix = np.clip(depth_matrix, None, clamp_threshold)
        
        # Convert to torch tensors
        self.depth = torch.tensor(depth_matrix, dtype=dtype, device=device)
        self.n_bins = self.depth.shape[0]
        self.n_samples = self.depth.shape[1]
        
        print(f"Loaded data: {self.n_bins} bins x {self.n_samples} samples")
        print(f"Depth range: [{self.depth.min():.3f}, {self.depth.max():.3f}]")
        print(f"Depth mean: {self.depth.mean():.3f}, median: {self.depth.median():.3f}")


class AneuploidyModel:
    """
    Hierarchical Bayesian model for aneuploidy detection from normalized read depth.
    
    Model structure:
    - Copy number (CN): discrete latent variable [0,1,2,3,4,5] with prior favoring CN=2
    - Per-bin mean bias: modulates expected depth at each bin
    - Per-sample variance: controls noise in each sample
    - Per-bin variance: controls noise at each bin
    - Observed depth ~ Normal(CN * sample_mean * bin_bias, variance)
    """
    
    def __init__(self,
                 n_states: int = 6,
                 alpha_ref: float = 50.0, 
                 alpha_non_ref: float = 1.0,
                 var_bias_bin: float = 0.1,
                 var_sample: float = 0.2,
                 var_bin: float = 0.2,
                 device: str = 'cpu',
                 dtype: torch.dtype = torch.float32,
                 debug: bool = False,
                 guide_type: str = 'diagonal'):
        """
        Args:
            n_states: Number of copy number states (default: 6 for [0,1,2,3,4,5])
            alpha_ref: Dirichlet concentration for reference state (CN=2)
            alpha_non_ref: Dirichlet concentration for non-reference states
            var_bias_bin: Variance for per-bin mean bias (log-normal)
            var_sample: Variance for per-sample variance factor (log-normal)
            var_bin: Variance for per-bin variance factor (log-normal)
            device: torch device
            dtype: torch data type
            debug: Whether to print debug statements in model()
            guide_type: Type of variational guide ('diagonal' for AutoDiagonalNormal, 'delta' for AutoDelta)
        """
        self.n_states = n_states
        self.alpha_ref = alpha_ref
        self.alpha_non_ref = alpha_non_ref
        self.var_bias_bin = var_bias_bin
        self.var_sample = var_sample
        self.var_bin = var_bin
        self.device = device
        self.dtype = dtype
        self.debug = debug
        self.guide_type = guide_type
        
        # Training history
        self.loss_history = {'epoch': [], 'elbo': []}
        
        # MCMC samples (if using MCMC)
        self.mcmc_samples = None
        self.inference_method = None  # 'svi' or 'mcmc'
        
        # Define which sites to expose to the guide (continuous latent variables)
        self.latent_sites = ['bin_bias', 'sample_var', 'bin_var', 'cn_probs']
        
        # Initialize guide based on type
        blocked_model = poutine.block(self.model, expose=self.latent_sites)
        if guide_type == 'delta':
            self.guide = AutoDelta(blocked_model)
        elif guide_type == 'diagonal':
            self.guide = AutoDiagonalNormal(blocked_model)
        else:
            raise ValueError(f"Unknown guide_type: {guide_type}. Choose 'diagonal' or 'delta'.")
    
    @config_enumerate(default="parallel")
    def model(self, depth: torch.Tensor, n_bins: int = None, n_samples: int = None):
        """
        Probabilistic model for aneuploidy detection.
        
        Args:
            depth: Observed normalized read depth (n_bins x n_samples)
        """
        
        if self.debug:
            print(f"\n=== MODEL DEBUG ===")
            print(f"depth.shape: {depth.shape}")
        
        zero_t = torch.zeros(1, device=self.device, dtype=self.dtype)
        one_t = torch.ones(1, device=self.device, dtype=self.dtype)
        cn_states = torch.arange(0, self.n_states, device=self.device, dtype=self.dtype)
        
        # Plates for bins and samples
        plate_bins = pyro.plate('bins', n_bins, dim=-2, device=self.device)
        plate_samples = pyro.plate('samples', n_samples, dim=-1, device=self.device)
        
        # Per-sample variance factor (log-normal prior)
        with plate_samples:
            sample_var = pyro.sample('sample_var', dist.Exponential(1.0 / self.var_sample))
        if self.debug:
            print(f"sample_var.shape: {sample_var.shape}")
        
        # Per-bin latent variables
        with plate_bins:
            # Per-bin mean bias factor (log-normal prior, centered at 1.0)
            bin_bias = pyro.sample('bin_bias', dist.LogNormal(zero_t, self.var_bias_bin))
            if self.debug:
                print(f"bin_bias.shape: {bin_bias.shape}")
            
            # Per-bin variance factor (log-normal prior)
            bin_var = pyro.sample('bin_var', dist.Exponential(1.0 / self.var_bin))
            if self.debug:
                print(f"bin_var.shape: {bin_var.shape}")
            
            # Copy number prior (Dirichlet-Categorical)
            # Heavily weight CN=2 (diploid)
            alpha_cn = self.alpha_non_ref * one_t.expand(self.n_states)
            alpha_cn[2] = self.alpha_ref  # Favor CN=2
            cn_probs = pyro.sample('cn_probs', dist.Dirichlet(alpha_cn))
        if self.debug:
            print(f"cn_probs.shape: {cn_probs.shape}")
        
        # Per-bin, per-sample copy number and observation
        with plate_bins, plate_samples:
            # Sample copy number (discrete latent variable)
            # Shape: (n_bins, n_samples)
            cn = pyro.sample('cn', dist.Categorical(cn_probs))
            if self.debug:
                print(f"cn.shape: {cn.shape}")
            
            # Expected depth: CN * bin_bias
            # For diploid (CN=2), expected depth should be ~2.0
            # bin_bias modulates this expectation per bin
            # CN=0 -> 0, CN=1 -> bin_bias, CN=2 -> 2*bin_bias, etc.
            if self.debug:
                print(f"bin_bias.shape: {bin_bias.shape}")
            expected_depth = Vindex(cn_states)[cn] * bin_bias
            if self.debug:
                print(f"expected_depth.shape: {expected_depth.shape}")
            
            # Variance: combination of sample and bin variance
            # Use summation to combine the two variance factors
            if self.debug:
                print(f"bin_var).shape: {bin_var.shape}")
            variance = sample_var + bin_var
            if self.debug:
                print(f"variance.shape: {variance.shape}")
            std = torch.sqrt(variance)
            if self.debug:
                print(f"std.shape: {std.shape}")
            
            # Observed depth
            if self.debug:
                print(f"About to sample obs with expected_depth.shape={expected_depth.shape}, std.shape={std.shape}, depth.shape={depth.shape}")
            pyro.sample('obs', dist.Normal(expected_depth, std), obs=depth)
        if self.debug:
            print(f"=== END MODEL DEBUG ===\n")
    
    def train(self,
              data,
              max_iter: int = 1000,
              lr_init: float = 0.01,
              lr_min: float = 0.001,
              lr_decay: float = 500,
              adam_beta1: float = 0.9,
              adam_beta2: float = 0.999,
              log_freq: int = 50,
              jit: bool = False):
        """
        Train the model using stochastic variational inference (SVI).
        
        Args:
            data: AneuploidyData object
            max_iter: Maximum number of training iterations
            lr_init: Initial learning rate
            lr_min: Minimum learning rate
            lr_decay: Learning rate decay constant
            adam_beta1: Adam optimizer beta1 parameter
            adam_beta2: Adam optimizer beta2 parameter
            log_freq: Frequency of logging (iterations)
            jit: Whether to use JIT compilation
        """
        print("Initializing training...")
        pyro.clear_param_store()
        
        # Create optimizer with learning rate schedule
        scheduler = pyro.optim.LambdaLR({
            'optimizer': torch.optim.Adam,
            'optim_args': {'lr': 1., 'betas': (adam_beta1, adam_beta2)},
            'lr_lambda': lambda k: lr_min + (lr_init - lr_min) * np.exp(-k / lr_decay)
        })
        
        # Create ELBO loss
        if jit:
            loss = JitTraceEnum_ELBO()
        else:
            loss = TraceEnum_ELBO()
        
        # Create SVI object
        svi = SVI(self.model, self.guide, optim=scheduler, loss=loss)
        
        print(f"Training for {max_iter} iterations...")
        print(f"Data: {data.n_bins} bins x {data.n_samples} samples")
        
        try:
            # Configure tqdm with mininterval to force updates
            with tqdm(range(max_iter), desc="Training", unit="epoch", 
                     mininterval=0.1, dynamic_ncols=True) as pbar:
                for epoch in pbar:
                    # Train step
                    epoch_loss = svi.step(depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples)
                    scheduler.step(epoch)
                    
                    # Record loss
                    self.loss_history['epoch'].append(epoch)
                    self.loss_history['elbo'].append(epoch_loss)
                    
                    # Update progress bar with loss
                    pbar.set_postfix({'loss': f'{epoch_loss:.4f}'})
                    
                    # Log progress
                    if (epoch + 1) % log_freq == 0:
                        tqdm.write(f"[epoch {epoch + 1:04d}]  loss: {epoch_loss:.4f}")
            
            print(f"\nTraining complete after {max_iter} epochs, final loss: {epoch_loss:.4f}")
            
        except KeyboardInterrupt:
            print("\nTraining interrupted by user.")
        
        self.inference_method = 'svi'
    
    def train_mcmc(self,
                   data,
                   num_samples: int = 1000,
                   warmup_steps: int = 200,
                   num_chains: int = 1,
                   step_size: float = 0.01,
                   kernel: str = 'NUTS'):
        """
        Train the model using MCMC sampling.
        
        Args:
            data: AneuploidyData object
            num_samples: Number of MCMC samples to draw
            warmup_steps: Number of warmup/burn-in steps
            num_chains: Number of MCMC chains
            step_size: Initial step size for HMC (ignored for NUTS)
            kernel: MCMC kernel to use ('NUTS' or 'HMC')
        """
        print(f"Initializing MCMC with {kernel} kernel...")
        print(f"  Warmup steps: {warmup_steps}")
        print(f"  Sampling steps: {num_samples}")
        print(f"  Chains: {num_chains}")
        pyro.clear_param_store()
        
        # Create a model without discrete sites for MCMC
        # We'll marginalize over CN using enumeration
        def continuous_model(depth, n_bins, n_samples):
            # This is a version of the model with CN marginalized out
            # For now, we'll sample it as part of MCMC
            return self.model(depth, n_bins, n_samples)
        
        # Choose kernel
        if kernel == 'NUTS':
            mcmc_kernel = NUTS(continuous_model, step_size=step_size, adapt_step_size=True,
                              target_accept_prob=0.8)
        elif kernel == 'HMC':
            mcmc_kernel = HMC(continuous_model, step_size=step_size, num_steps=10)
        else:
            raise ValueError(f"Unknown kernel: {kernel}. Choose 'NUTS' or 'HMC'.")
        
        # Create MCMC object
        mcmc = MCMC(mcmc_kernel, num_samples=num_samples, warmup_steps=warmup_steps,
                    num_chains=num_chains, disable_progbar=False)
        
        # Run MCMC
        print("Running MCMC sampling...")
        mcmc.run(depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples)
        
        # Store samples
        self.mcmc_samples = mcmc.get_samples()
        self.inference_method = 'mcmc'
        
        # Print diagnostics
        print("\nMCMC Diagnostics:")
        mcmc.summary(prob=0.95)
        
        print("\nMCMC sampling complete.")
    
    def get_map_estimates(self, data):
        """
        Get MAP (maximum a posteriori) estimates of all latent variables.
        
        Returns:
            Dictionary with MAP estimates
        """
        print("Computing MAP estimates...")
        
        if self.inference_method == 'mcmc':
            # For MCMC, use posterior mean as point estimate
            print("Using posterior mean from MCMC samples...")
            map_estimates = {}
            
            for key in ['bin_bias', 'sample_var', 'bin_var', 'cn_probs']:
                if key in self.mcmc_samples:
                    samples = self.mcmc_samples[key]
                    map_estimates[key] = samples.mean(axis=0).detach().cpu().numpy()
            
            # For CN, use mode (most common value) across samples
            if 'cn' in self.mcmc_samples:
                cn_samples = self.mcmc_samples['cn'].detach().cpu().numpy()
                # Compute mode along sample dimension (axis=0)
                from scipy import stats
                map_estimates['cn'] = stats.mode(cn_samples, axis=0, keepdims=False)[0]
            
        else:
            # SVI mode: use guide parameters
            # Get guide trace (contains MAP estimates of continuous variables)
            guide_trace = poutine.trace(self.guide).get_trace(depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples)
            
            # Get model trace conditioned on guide
            trained_model = poutine.replay(self.model, trace=guide_trace)
            
            # Get discrete MAP using infer_discrete
            inferred_model = infer_discrete(trained_model, temperature=0, first_available_dim=-3)
            trace = poutine.trace(inferred_model).get_trace(depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples)
            
            # Extract all latent variables
            map_estimates = {}
            
            # Continuous variables from guide
            map_estimates['bin_bias'] = guide_trace.nodes['bin_bias']['value'].detach().cpu().numpy()
            map_estimates['sample_var'] = guide_trace.nodes['sample_var']['value'].detach().cpu().numpy()
            map_estimates['bin_var'] = guide_trace.nodes['bin_var']['value'].detach().cpu().numpy()
            map_estimates['cn_probs'] = guide_trace.nodes['cn_probs']['value'].detach().cpu().numpy()
            
            # Discrete variable (copy number)
            map_estimates['cn'] = trace.nodes['cn']['value'].detach().cpu().numpy()
        
        print("MAP estimates computed.")
        return map_estimates
    
    def run_discrete_inference(self, data, n_samples: int = 1000, log_freq: int = 100):
        """
        Run discrete inference to get posterior distribution over copy numbers.
        
        Args:
            data: AneuploidyData object
            n_samples: Number of samples for discrete inference
            log_freq: Logging frequency
            
        Returns:
            Dictionary with copy number posterior probabilities
        """
        print(f"Running discrete inference with {n_samples} samples...")
        
        # Get guide trace
        guide_trace = poutine.trace(self.guide).get_trace(depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples)
        trained_model = poutine.replay(self.model, trace=guide_trace)
        
        # Accumulate samples
        cn_samples = []
        
        with torch.no_grad():
            for i in range(n_samples):
                inferred_model = infer_discrete(trained_model, temperature=1, first_available_dim=-3)
                trace = poutine.trace(inferred_model).get_trace(depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples)
                cn = trace.nodes['cn']['value'].detach().cpu().numpy()
                cn_samples.append(cn)
                
                if (i + 1) % log_freq == 0:
                    print(f"[sample {i + 1:04d}]")
        
        # Stack samples and compute frequencies
        cn_samples = np.array(cn_samples)  # (n_samples, n_bins, n_samples)
        
        # Convert to one-hot encoding and compute frequencies
        cn_freq = np.zeros((data.n_bins, data.n_samples, self.n_states))
        for state in range(self.n_states):
            cn_freq[..., state] = (cn_samples == state).mean(axis=0)
        
        print("Discrete inference complete.")
        return {'cn_posterior': cn_freq}
    

def plot_training_loss(model: AneuploidyModel, output_dir: str):
    """Plot training loss over epochs"""
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(model.loss_history['epoch'], model.loss_history['elbo'])
    ax.set_xlabel('Epoch')
    ax.set_ylabel('ELBO')
    ax.set_title('Training Loss')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    output_path = os.path.join(output_dir, 'training_loss.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {output_path}")
    plt.close()


def add_chromosome_labels(ax, chr_array, x_transformed=None):
    """
    Add chromosome labels to x-axis of a plot.
    
    Args:
        ax: matplotlib axis
        chr_array: array of chromosome names for each bin
        x_transformed: optional array of transformed x coordinates (for equal-width chromosomes)
    """
    # Find chromosome boundaries
    chr_changes = np.where(chr_array[:-1] != chr_array[1:])[0] + 1
    chr_boundaries = np.concatenate([[0], chr_changes, [len(chr_array)]])
    
    # Calculate midpoint of each chromosome for labels
    chr_labels = []
    chr_positions = []
    chr_boundary_positions = []
    
    for i in range(len(chr_boundaries) - 1):
        start = chr_boundaries[i]
        end = chr_boundaries[i + 1]
        
        if x_transformed is not None:
            # Use transformed coordinates - chromosome centers are at i+0.5
            chr_positions.append(i + 0.5)
            if i > 0:
                chr_boundary_positions.append(i)
        else:
            # Use original bin indices
            mid = (start + end) / 2
            chr_positions.append(mid)
            if start > 0:
                chr_boundary_positions.append(start)
        
        chr_name = chr_array[start]
        chr_labels.append(str(chr_name).replace('chr', ''))
    
    # Add vertical lines at chromosome boundaries
    for boundary in chr_boundary_positions:
        ax.axvline(boundary, color='gray', linestyle='--', alpha=1, linewidth=1)
    
    # Set x-axis labels
    ax.set_xticks(chr_positions)
    ax.set_xticklabels(chr_labels, rotation=0, ha='center')
    ax.set_xlabel('Chromosome')


def plot_sample(data: AneuploidyData, map_estimates: dict, cn_posterior: dict, output_dir: str,
                sample_idx: int = 0, is_aneuploid: bool = False, aneuploid_chrs: list = None):
    """
    Combined plot showing MAP estimates and CN posterior.
    
    Args:
        data: AneuploidyData object
        map_estimates: Dictionary with MAP estimates
        cn_posterior: Dictionary with 'cn_posterior' key
        sample_idx: Index of sample to plot
        is_aneuploid: Whether this sample has high confidence aneuploidy
        aneuploid_chrs: List of tuples (chr_name, cn, prob) for aneuploid chromosomes
    """
    fig, axes = plt.subplots(3, 1, figsize=(8, 8))
    
    # Add sample name as figure title positioned at left edge
    fig.text(0.1/8, 0.98, f'{data.sample_ids[sample_idx]}', 
             fontsize=12, fontweight='bold', va='top', ha='left')
    
    # Extract data
    observed = data.depth[:, sample_idx].cpu().numpy()
    cn = map_estimates['cn'][:, sample_idx]
    sample_var_all = map_estimates['sample_var'].flatten()  # All samples' variance factors
    sample_var_current = sample_var_all[sample_idx]  # Current sample's variance factor
    cn_probs = cn_posterior['cn_posterior'][:, sample_idx, :]  # (n_bins, n_states)
    
    colors = [ '#004D40', '#FFC107',  '#1E88E5', '#D81B60', '#38006B']
    
    # Create normalized x-axis where each chromosome has equal width
    # Get chromosome boundaries
    chr_array = data.chr
    chr_changes = np.where(chr_array[:-1] != chr_array[1:])[0] + 1
    chr_boundaries = np.concatenate([[0], chr_changes, [len(chr_array)]])
    
    # Calculate transformed x positions (each chromosome gets unit width)
    x_transformed = np.zeros(len(observed))
    for i in range(len(chr_boundaries) - 1):
        start_idx = chr_boundaries[i]
        end_idx = chr_boundaries[i + 1]
        n_bins_in_chr = end_idx - start_idx
        # Map bins within this chromosome to interval [i, i+1]
        x_transformed[start_idx:end_idx] = i + np.linspace(0, 1, n_bins_in_chr, endpoint=False)
    
    # Plot 1: Observed depth with MAP Copy number calls overlaid
    ax = axes[0]
    # First plot MAP copy number calls
    ax.plot(x_transformed, cn, 'o', alpha=0.5, markersize=4,
            label=f'Predicted bin copy number', color='black')
    # Then plot observed depth on top so it's visible
    ax.plot(x_transformed, observed, 'r-', linewidth=1, alpha=0.7, label='Normalized read depth')
    ax.legend(loc='lower right', bbox_to_anchor=(1.0, 1.02), ncol=2, borderaxespad=0)
    ax.grid(True, axis='y', alpha=1, linestyle='-', linewidth=1)
    ax.set_ylim([-0.5, 4.5])
    ax.set_xlim([x_transformed.min(), x_transformed.max()])
    
    # Plot 2: Stacked area plot of CN probabilities
    ax = axes[1]
    ax.stackplot(x_transformed, cn_probs.T, labels=[f'CN={i}' for i in range(6)], 
                alpha=0.7, colors=colors)
    ax.set_ylabel('Copy Number Probability')
    ax.legend(loc='lower right', bbox_to_anchor=(1.0, 1.02), ncol=6, borderaxespad=0)
    ax.set_ylim([0, 1])
    ax.set_xlim([x_transformed.min(), x_transformed.max()])
    
    # Plot 3: Sample variance distribution
    ax = axes[2]
    ax.hist(np.sqrt(sample_var_all), bins=30, alpha=0.7, edgecolor='black', linewidth=0.5, color='gray')
    ax.axvline(np.sqrt(sample_var_current), color='red', linestyle='--', linewidth=2, 
              label=f'This sample: {np.sqrt(sample_var_current):.3f}')
    ax.set_xlabel('Sample stdev')
    ax.set_ylabel('Count')
    ax.set_yscale('log')
    ax.legend(loc='lower right', bbox_to_anchor=(1.0, 1.02), borderaxespad=0)
    ax.grid(True, alpha=0.3)
    
    # Highlight aneuploid chromosomes in first 3 panels if applicable
    if is_aneuploid and aneuploid_chrs:
        # Extract chromosome names from aneuploid_chrs
        aneuploid_chr_names = set([chr_name for chr_name, _, _ in aneuploid_chrs])
        
        # Find bin ranges for each aneuploid chromosome
        for chr_name in aneuploid_chr_names:
            chr_mask = data.chr == chr_name
            chr_indices = np.where(chr_mask)[0]
            
            if len(chr_indices) > 0:
                chr_start_idx = chr_indices[0]
                chr_end_idx = chr_indices[-1]
                
                # Get transformed coordinates for highlighting
                chr_start_transformed = x_transformed[chr_start_idx]
                chr_end_transformed = x_transformed[chr_end_idx]
                
                # Add red transparent box to first 2 panels
                for ax in axes[:2]:
                    y_min, y_max = ax.get_ylim()
                    ax.axvspan(chr_start_transformed, chr_end_transformed, alpha=0.15, color='red', zorder=0)
    
    # Add chromosome labels to first 2 panels
    for ax in axes[:2]:
        add_chromosome_labels(ax, data.chr, x_transformed=x_transformed)
    
    plt.tight_layout()
    sample_name = data.sample_ids[sample_idx].replace('/', '_').replace(' ', '_')
    aneu_suffix = 'ANEU' if is_aneuploid else 'NORMAL'
    filename = f'sample_{aneu_suffix}_{sample_name}.png'
    
    # Create subdirectory based on aneuploidy status
    subdir = 'sample_plots'
    subdir_path = os.path.join(output_dir, subdir)
    os.makedirs(subdir_path, exist_ok=True)
    
    output_path = os.path.join(subdir_path, filename)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {output_path}")
    plt.close()


def plot_bin_variance_bias(data: AneuploidyData, map_estimates: dict, output_dir: str):
    """
    Plot bin bias and bin variance posteriors (global across all samples).
    
    Args:
        data: AneuploidyData object
        map_estimates: Dictionary with MAP estimates
    """
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))
    
    # Extract data
    bin_bias = map_estimates['bin_bias'].flatten()
    bin_var = map_estimates['bin_var'].flatten()
    
    x = np.arange(len(bin_bias))
    
    # Plot 1: Bin bias
    ax = axes[0]
    ax.plot(x, bin_bias, 'o', alpha=0.5, markersize=3, color='black', label='Bin bias')
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='No bias')
    ax.set_ylabel('Bin mean bias')
    ax.set_title('Bin Posteriors')
    ax.set_xlim([x.min(), x.max()])
    
    # Plot 2: Bin variance
    ax = axes[1]
    ax.plot(x, np.sqrt(bin_var), 'o', alpha=0.5, markersize=3, color='purple', label='Bin stdev')
    ax.set_ylabel('Bin stdev')
    ax.set_xlim([x.min(), x.max()])
    
    # Add chromosome labels to both panels
    for ax in axes:
        add_chromosome_labels(ax, data.chr)
    
    plt.tight_layout()
    filename = 'bin_posteriors.png'
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {output_path}")
    plt.close()


def plot_map_estimates(data: AneuploidyData, map_estimates: dict, sample_idx: int = 0):
    """
    Plot MAP estimates for a single sample.
    
    Args:
        data: AneuploidyData object
        map_estimates: Dictionary with MAP estimates
        sample_idx: Index of sample to plot
    """
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    
    # Extract data for this sample
    observed = data.depth[:, sample_idx].cpu().numpy()
    cn = map_estimates['cn'][:, sample_idx]
    bin_bias = map_estimates['bin_bias']
    
    # Plot 1: Observed depth and copy number
    ax = axes[0]
    ax.plot(range(len(observed)), observed, 'o-', alpha=0.5, markersize=3, linewidth=0.5, label='Observed depth')
    ax.set_ylabel('Normalized depth')
    ax.set_title(f'Sample: {data.sample_ids[sample_idx]}')
    ax.legend()
    ax.grid(True, alpha=0.5)
    
    # Plot 2: Copy number calls
    ax = axes[1]
    colors = ['red', 'orange', 'green', 'blue', 'purple', 'brown']
    for state in range(6):
        mask = cn == state
        if np.any(mask):
            indices = np.where(mask)[0]
            ax.plot(indices, cn[mask], 'o-', alpha=0.7, markersize=4, linewidth=1,
                   label=f'CN={state}', color=colors[state])
    ax.set_ylabel('Copy Number')
    ax.set_ylim(-0.5, 5.5)
    ax.legend()
    ax.grid(True, alpha=0.5)
    
    # Plot 3: Bin bias
    ax = axes[2]
    ax.plot(range(len(bin_bias)), bin_bias, 'o-', alpha=0.5, markersize=3, linewidth=0.5, color='black', label='Bin bias')
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='No bias')
    ax.set_ylabel('Bin bias')
    ax.legend()
    ax.grid(True, alpha=0.5)
    
    # Add chromosome labels to all panels
    for ax in axes:
        add_chromosome_labels(ax, data.chr)
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, f'map_estimates_sample_{sample_idx}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {output_path}")
    plt.close()


def plot_cn_posterior(data: AneuploidyData, cn_posterior: dict, sample_idx: int = 0):
    """
    Plot copy number posterior probabilities and observed depth.
    
    Args:
        data: AneuploidyData object
        cn_posterior: Dictionary with 'cn_posterior' key
        sample_idx: Index of sample to plot
    """
    cn_probs = cn_posterior['cn_posterior'][:, sample_idx, :]  # (n_bins, n_states)
    
    # Get observed depth
    obs_actual = data.depth[:, sample_idx].cpu().numpy()  # (n_bins,)
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    
    # Plot 1: Stacked area plot of CN probabilities
    ax = axes[0]
    colors = ['red', 'orange', 'green', 'blue', 'purple', 'brown']
    x = np.arange(cn_probs.shape[0])
    
    ax.stackplot(x, cn_probs.T, labels=[f'CN={i}' for i in range(6)], 
                alpha=0.7, colors=colors)
    ax.set_ylabel('Probability')
    ax.set_title(f'Copy Number Posterior - Sample: {data.sample_ids[sample_idx]}')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Observed depth
    ax = axes[1]
    ax.plot(x, obs_actual, 'k-', linewidth=0.8, alpha=0.7, label='Observed')
    ax.set_ylabel('Normalized read depth')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add chromosome labels to both panels
    for ax in axes:
        add_chromosome_labels(ax, data.chr)
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, f'cn_posterior_sample_{sample_idx}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {output_path}")
    plt.close()


def summary_statistics(data: AneuploidyData, map_estimates: dict, cn_posterior: dict, 
                      min_fraction: float = 0.5, prob_threshold: float = 0.5):
    """
    Print summary statistics of the inference results.
    
    Args:
        data: AneuploidyData object
        map_estimates: Dictionary with MAP estimates
        cn_posterior: Dictionary with copy number posterior
        min_fraction: Minimum fraction of bins that must be abnormal to call chromosome aneuploid
        prob_threshold: Probability threshold for high-confidence calls
    """
    print("="*80)
    print("INFERENCE SUMMARY")
    print("="*80)
    
    # Sample variance statistics
    sample_var = map_estimates['sample_var']
    print(f"\nSample variance factors:")
    print(f"  Mean: {sample_var.mean():.4f}")
    print(f"  Std:  {sample_var.std():.4f}")
    print(f"  Min:  {sample_var.min():.4f}")
    print(f"  Max:  {sample_var.max():.4f}")
    
    # Bin variance statistics
    bin_var = map_estimates['bin_var']
    print(f"\nBin variance factors:")
    print(f"  Mean: {bin_var.mean():.4f}")
    print(f"  Std:  {bin_var.std():.4f}")
    print(f"  Min:  {bin_var.min():.4f}")
    print(f"  Max:  {bin_var.max():.4f}")
    
    # Bin bias statistics
    bin_bias = map_estimates['bin_bias']
    print(f"\nBin bias factors:")
    print(f"  Mean: {bin_bias.mean():.4f}")
    print(f"  Std:  {bin_bias.std():.4f}")
    print(f"  Min:  {bin_bias.min():.4f}")
    print(f"  Max:  {bin_bias.max():.4f}")
    
    # Copy number statistics
    cn = map_estimates['cn']
    cn_probs = cn_posterior['cn_posterior']
    
    print(f"\nCopy number calls (MAP):")
    for state in range(6):
        count = (cn == state).sum()
        pct = 100 * count / cn.size
        print(f"  CN={state}: {count:6d} calls ({pct:5.2f}%)")
    
    # Per-chromosome aneuploidy detection
    print(f"\n" + "="*80)
    print(f"PER-CHROMOSOME ANEUPLOIDY DETECTION")
    print(f"  Probability threshold: {prob_threshold}")
    print("="*80)
    
    cn_probs = cn_posterior['cn_posterior']
    
    # Get unique chromosomes
    unique_chrs = np.unique(data.chr)
    
    # Separate sex chromosomes from autosomes
    sex_chrs = ['chrX', 'chrY']
    autosomes = [chr_name for chr_name in unique_chrs if chr_name not in sex_chrs]
    
    # Track aneuploid chromosomes per sample
    aneuploid_samples = {}
    for sample_idx in range(data.n_samples):
        aneuploid_samples[sample_idx] = []
    
    # Helper function to get chromosome CN and mean probability
    def get_chr_info(chr_name, sample_idx):
        chr_mask = data.chr == chr_name
        chr_bins = np.sum(chr_mask)
        if chr_bins == 0:
            return None, 0, 0
        
        # Determine most common CN state
        chr_cn = cn[chr_mask, sample_idx]
        cn_counts = np.bincount(chr_cn, minlength=6)
        most_common_cn = np.argmax(cn_counts)
        
        # Calculate mean probability of the called CN state across all bins
        chr_cn_probs = cn_probs[chr_mask, sample_idx, most_common_cn]
        mean_prob = chr_cn_probs.mean()
        
        return most_common_cn, mean_prob, chr_bins
    
    print(f"\nChromosome-level aneuploidy calls:")
    
    # Process autosomes
    for chr_name in autosomes:
        chr_mask = data.chr == chr_name
        chr_bins = np.sum(chr_mask)
        
        print(f"\n  {chr_name} ({chr_bins} bins):")
        
        # Check each sample for this chromosome
        for sample_idx in range(data.n_samples):
            sample_name = data.sample_ids[sample_idx]
            
            # Get CN and mean probability for this chromosome
            chr_cn_state, chr_mean_prob, _ = get_chr_info(chr_name, sample_idx)
            
            # For autosomes, call aneuploidy if CN != 2 and has high confidence
            is_aneuploid = (chr_cn_state != 2) and (chr_mean_prob > prob_threshold)
            
            if is_aneuploid:
                print(f"    {sample_name}: ANEUPLOID (CN≈{chr_cn_state}, P(CN={chr_cn_state})={chr_mean_prob:.2f})")
                aneuploid_samples[sample_idx].append((chr_name, chr_cn_state, chr_mean_prob))
    
    # Process sex chromosomes jointly
    print(f"\n  Sex chromosomes (chrX and chrY):")
    for sample_idx in range(data.n_samples):
        sample_name = data.sample_ids[sample_idx]
        
        # Get CN and mean probability for both sex chromosomes
        chrX_cn, chrX_mean_prob, chrX_bins = get_chr_info('chrX', sample_idx)
        chrY_cn, chrY_mean_prob, chrY_bins = get_chr_info('chrY', sample_idx)
        
        # Skip if both chromosomes are missing
        if chrX_cn is None and chrY_cn is None:
            continue
        
        # Handle case where one chromosome is missing
        if chrX_cn is None:
            chrX_cn, chrX_mean_prob, chrX_bins = 0, 0, 0
        if chrY_cn is None:
            chrY_cn, chrY_mean_prob, chrY_bins = 0, 0, 0
        
        # For sex chromosomes, check if we have high confidence in the CN call
        chrX_high_conf = chrX_mean_prob > prob_threshold if chrX_bins > 0 else True
        chrY_high_conf = chrY_mean_prob > prob_threshold if chrY_bins > 0 else True
        
        # Determine karyotype
        # Normal: XX (chrX=2, chrY=0) or XY (chrX=1, chrY=1)
        is_XX = (chrX_cn == 2 and chrY_cn == 0)
        is_XY = (chrX_cn == 1 and chrY_cn == 1)
        is_normal = is_XX or is_XY
        
        # Call sex chromosome aneuploidy if karyotype is abnormal and both have high confidence
        if not is_normal and chrX_high_conf and chrY_high_conf:
            karyotype = f"chrX={chrX_cn}, chrY={chrY_cn}"
            prob_str = f"P(chrX={chrX_cn})={chrX_mean_prob:.2f}, P(chrY={chrY_cn})={chrY_mean_prob:.2f}"
            print(f"    {sample_name}: SEX CHROMOSOME ANEUPLOIDY ({karyotype}, {prob_str})")
            
            # Add both chromosomes to the aneuploid list
            if chrX_bins > 0:
                aneuploid_samples[sample_idx].append(('chrX', chrX_cn, chrX_mean_prob))
            if chrY_bins > 0:
                aneuploid_samples[sample_idx].append(('chrY', chrY_cn, chrY_mean_prob))
        elif not is_normal:
            # One or both chromosomes don't meet threshold - report but don't call aneuploidy
            karyotype = f"chrX={chrX_cn} (P={chrX_mean_prob:.2f}), chrY={chrY_cn} (P={chrY_mean_prob:.2f})"
            print(f"    {sample_name}: Uncertain ({karyotype})")
    
    # Summary per sample
    print(f"\n" + "="*80)
    print(f"ANEUPLOIDY SUMMARY BY SAMPLE")
    print("="*80)
    for sample_idx in range(data.n_samples):
        sample_name = data.sample_ids[sample_idx]
        aneuploid_chrs = aneuploid_samples[sample_idx]
        if len(aneuploid_chrs) > 1:
            print(f"\n{sample_name}: ANEUPLOID ({len(aneuploid_chrs)} chromosome(s))")
            for chr_name, cn, mean_prob in aneuploid_chrs:
                print(f"  - {chr_name}: CN≈{cn} (P(CN={cn})={mean_prob:.2f})")
    
    print("="*80)
    
    return aneuploid_samples


def read_data(file_path: str) -> pd.DataFrame:
    
    # Load file
    print(f"Loading: {file_path}")
    try:
        # Read gzipped TSV file
        df_file = pd.read_csv(file_path, sep='\t', compression='infer')
        
        # Get source file identifier
        source_file = str(Path(file_path).parent.parent.name)
        
        # Rename chromosome column if it exists
        if '#Chr' in df_file.columns:
            df_file['Chr'] = df_file['#Chr']
            df_file = df_file.drop('#Chr', axis=1)
        
        # Create bin identifier
        df_file['Bin'] = df_file['Chr'].astype(str) + ':' + df_file['Start'].astype(str) + '-' + df_file['End'].astype(str)
        df_file = df_file.set_index('Bin')
        df_file['source_file'] = source_file
        return df_file
    except Exception as e:
        print(f"Error loading {file_path}")
        raise e


def filter_low_quality_bins(df: pd.DataFrame,
                           autosome_median_min: float = 1.5,
                           autosome_median_max: float = 2.5,
                           autosome_mad_max: float = 0.5,
                           chrX_median_min: float = 0.75,
                           chrX_median_max: float = 2.5,
                           chrX_mad_max: float = 0.5,
                           chrY_median_min: float = 0.0,
                           chrY_median_max: float = 2.0,
                           chrY_mad_max: float = 0.5):
    """
    Filter out low quality bins based on median and MAD thresholds.
    
    Args:
        df: DataFrame with bins as rows and samples as columns
        autosome_median_min: Minimum median depth for autosomal bins
        autosome_median_max: Maximum median depth for autosomal bins
        autosome_mad_max: Maximum MAD for autosomal bins
        chrX_median_min: Minimum median depth for chrX bins
        chrX_median_max: Maximum median depth for chrX bins
        chrX_mad_max: Maximum MAD for chrX bins
        chrY_median_min: Minimum median depth for chrY bins
        chrY_median_max: Maximum median depth for chrY bins
        chrY_mad_max: Maximum MAD for chrY bins
        
    Returns:
        Filtered DataFrame
    """
    # Get sample columns (exclude metadata)
    metadata_cols = ['Chr', 'Start', 'End', 'source_file']
    sample_cols = [col for col in df.columns if col not in metadata_cols]
    
    # Compute median and MAD for each bin across samples
    depths = df[sample_cols].values
    medians = np.median(depths, axis=1)
    mads = np.median(np.abs(depths - medians[:, np.newaxis]), axis=1)
    
    print(f"\n{'='*80}")
    print("FILTERING LOW QUALITY BINS")
    print(f"{'='*80}")
    print(f"Starting bins: {len(df)}")
    
    # Create filter masks for each chromosome type
    autosome_mask = ~df['Chr'].isin(['chrX', 'chrY'])
    chrX_mask = df['Chr'] == 'chrX'
    chrY_mask = df['Chr'] == 'chrY'
    
    # Initialize keep mask (all True)
    keep_mask = np.ones(len(df), dtype=bool)
    
    # Filter autosomes
    if autosome_mask.any():
        autosome_keep = (
            (medians >= autosome_median_min) &
            (medians <= autosome_median_max) &
            (mads <= autosome_mad_max)
        )
        autosome_filtered = autosome_mask & ~autosome_keep
        keep_mask[autosome_mask] = autosome_keep[autosome_mask]
        print(f"\nAutosomes:")
        print(f"  Thresholds: median [{autosome_median_min}, {autosome_median_max}], MAD <= {autosome_mad_max}")
        print(f"  Bins before: {autosome_mask.sum()}")
        print(f"  Bins filtered: {autosome_filtered.sum()}")
        print(f"  Bins after: {(autosome_mask & keep_mask).sum()}")
    
    # Filter chrX
    if chrX_mask.any():
        chrX_keep = (
            (medians >= chrX_median_min) &
            (medians <= chrX_median_max) &
            (mads <= chrX_mad_max)
        )
        chrX_filtered = chrX_mask & ~chrX_keep
        keep_mask[chrX_mask] = chrX_keep[chrX_mask]
        print(f"\nchrX:")
        print(f"  Thresholds: median [{chrX_median_min}, {chrX_median_max}], MAD <= {chrX_mad_max}")
        print(f"  Bins before: {chrX_mask.sum()}")
        print(f"  Bins filtered: {chrX_filtered.sum()}")
        print(f"  Bins after: {(chrX_mask & keep_mask).sum()}")
    
    # Filter chrY
    if chrY_mask.any():
        chrY_keep = (
            (medians >= chrY_median_min) &
            (medians <= chrY_median_max) &
            (mads <= chrY_mad_max)
        )
        chrY_filtered = chrY_mask & ~chrY_keep
        keep_mask[chrY_mask] = chrY_keep[chrY_mask]
        print(f"\nchrY:")
        print(f"  Thresholds: median [{chrY_median_min}, {chrY_median_max}], MAD <= {chrY_mad_max}")
        print(f"  Bins before: {chrY_mask.sum()}")
        print(f"  Bins filtered: {chrY_filtered.sum()}")
        print(f"  Bins after: {(chrY_mask & keep_mask).sum()}")
    
    # Apply filter
    df_filtered = df[keep_mask].copy()
    
    print(f"\nTotal bins after filtering: {len(df_filtered)}")
    print(f"Total bins removed: {len(df) - len(df_filtered)}")
    print(f"{'='*80}\n")
    
    return df_filtered


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Bayesian aneuploidy detection from normalized read depth data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input/Output
    parser.add_argument('-i', '--input', required=True,
                        help='Input TSV file with normalized read depth (bins x samples)')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Output directory for results and plots')
    
    # Model parameters
    parser.add_argument('--alpha-ref', type=float, default=1.0,
                        help='Dirichlet concentration parameter for reference copy number state (CN=2)')
    parser.add_argument('--alpha-non-ref', type=float, default=1.0,
                        help='Dirichlet concentration parameter for non-reference copy number states')
    parser.add_argument('--var-bias-bin', type=float, default=0.01,
                        help='Variance parameter for per-bin mean bias (log-normal prior)')
    parser.add_argument('--var-sample', type=float, default=0.01,
                        help='Variance parameter for per-sample variance factor')
    parser.add_argument('--var-bin', type=float, default=0.01,
                        help='Variance parameter for per-bin variance factor')
    parser.add_argument('--guide-type', type=str, default='delta', choices=['delta', 'diagonal'],
                        help='Type of variational guide (delta=MAP inference, diagonal=mean-field)')
    
    # Training parameters
    parser.add_argument('--max-iter', type=int, default=1000,
                        help='Maximum number of training iterations')
    parser.add_argument('--lr-init', type=float, default=0.02,
                        help='Initial learning rate')
    parser.add_argument('--lr-min', type=float, default=0.01,
                        help='Minimum learning rate')
    parser.add_argument('--lr-decay', type=float, default=500,
                        help='Learning rate decay constant')
    parser.add_argument('--log-freq', type=int, default=50,
                        help='Logging frequency (iterations)')
    parser.add_argument('--jit', action='store_true', default=False,
                        help='Enable JIT compilation for faster training')
    
    # Inference parameters
    parser.add_argument('--n-discrete-samples', type=int, default=1000,
                        help='Number of samples for discrete inference')
    parser.add_argument('--prob-threshold', type=float, default=0.5,
                        help='Probability threshold for high-confidence aneuploidy calls')
    
    # Data filtering
    parser.add_argument('--viable-only', action='store_true', default=False,
                        help='Subset data to viable trisomy chromosomes only (chr13, chr18, chr21, chrX, chrY)')
    
    # Bin quality filtering parameters
    parser.add_argument('--autosome-median-min', type=float, default=1.75,
                        help='Minimum median depth for autosomal bins')
    parser.add_argument('--autosome-median-max', type=float, default=2.25,
                        help='Maximum median depth for autosomal bins')
    parser.add_argument('--autosome-mad-max', type=float, default=0.25,
                        help='Maximum MAD for autosomal bins')
    parser.add_argument('--chrX-median-min', type=float, default=0.75,
                        help='Minimum median depth for chrX bins')
    parser.add_argument('--chrX-median-max', type=float, default=2.25,
                        help='Maximum median depth for chrX bins')
    parser.add_argument('--chrX-mad-max', type=float, default=0.5,
                        help='Maximum MAD for chrX bins')
    parser.add_argument('--chrY-median-min', type=float, default=0.0,
                        help='Minimum median depth for chrY bins')
    parser.add_argument('--chrY-median-max', type=float, default=1.25,
                        help='Maximum median depth for chrY bins')
    parser.add_argument('--chrY-mad-max', type=float, default=0.5,
                        help='Maximum MAD for chrY bins')
    
    # Device
    parser.add_argument('--device', type=str, default='cpu', choices=['cpu', 'cuda'],
                        help='Device to use for computation')
    
    return parser.parse_args()


def main():
    """
    Main function to run aneuploidy detection pipeline.
    """
    # Parse command-line arguments
    args = parse_args()
    
    bins_pattern = args.input
    output_dir = args.output_dir

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Delete all existing sample_*.png files
    sample_pattern = os.path.join(output_dir, "sample_plots/*.png")
    existing_files = glob.glob(sample_pattern)
    if existing_files:
        print(f"Deleting {len(existing_files)} existing sample_plots/*.png files...")
        for file_path in existing_files:
            os.remove(file_path)
        print("Cleanup complete.")

    df = read_data(bins_pattern)
    
    # Filter to viable trisomy chromosomes if requested
    if args.viable_only:
        viable_chrs = {'chr13', 'chr18', 'chr21', 'chrX', 'chrY'}
        n_bins_before = len(df)
        df = df[df['Chr'].isin(viable_chrs)]
        n_bins_after = len(df)
        print(f"Filtered to viable chromosomes only: {n_bins_before} -> {n_bins_after} bins")
        print(f"Chromosomes retained: {sorted(df['Chr'].unique())}")

    # Filter low quality bins
    df = filter_low_quality_bins(df,
                                autosome_median_min=args.autosome_median_min,
                                autosome_median_max=args.autosome_median_max,
                                autosome_mad_max=args.autosome_mad_max,
                                chrX_median_min=args.chrX_median_min,
                                chrX_median_max=args.chrX_median_max,
                                chrX_mad_max=args.chrX_mad_max,
                                chrY_median_min=args.chrY_median_min,
                                chrY_median_max=args.chrY_median_max,
                                chrY_mad_max=args.chrY_mad_max)

    print(df)

    # Set up Pyro
    pyro.enable_validation(True)
    pyro.distributions.enable_validation(True)
    pyro.set_rng_seed(42)
    torch.manual_seed(42)
    np.random.seed(42)


    # Prepare data - optionally subset for faster testing
    # Uncomment the line below to test on chr X only
    # df_subset = df[df['Chr'] == 'chrX'].copy()
    # Or use all data:
    df_subset = df.copy()

    # Create data object
    device = args.device  # Use 'cuda' if GPU available
    dtype = torch.float32
    inference_method = 'svi'

    data = AneuploidyData(df_subset, device=device, dtype=dtype)

    # Initialize model
    model = AneuploidyModel(
        n_states=4,                        # CN states (hard-coded to 0,1,2,3)
        alpha_ref=args.alpha_ref,         # Prior on CN=2
        alpha_non_ref=args.alpha_non_ref,     # Prior on other states
        var_bias_bin=args.var_bias_bin,     # Bin-to-bin bias variation
        var_sample=args.var_sample,       # Sample variance
        var_bin=args.var_bin,          # Bin variance
        device=device,
        dtype=dtype,
        guide_type=args.guide_type      # Variational guide type
    )

    # Train model
    if inference_method == 'svi':
        model.train(
            data=data,
            max_iter=args.max_iter,         # Number of training iterations
            lr_init=args.lr_init,         # Initial learning rate
            lr_min=args.lr_min,         # Minimum learning rate
            lr_decay=args.lr_decay,         # Learning rate decay constant
            log_freq=args.log_freq,          # Log frequency
            jit=args.jit             # JIT compilation
        )
    else:
        model.train_mcmc(
            data=data,
            warmup_steps=100,        # Number of warmup iterations
            num_samples=200,       # Number of posterior samples to draw
            num_chains=1,           # Number of MCMC chains
            step_size=0.01,
            kernel='NUTS'
        )

    # Plot training loss
    plot_training_loss(model, output_dir=output_dir)

    # Get MAP estimates
    map_estimates = model.get_map_estimates(data)

    # Run discrete inference to get CN posterior probabilities
    cn_posterior = model.run_discrete_inference(data, n_samples=args.n_discrete_samples, log_freq=args.log_freq)
    print("cn_posterior keys:", cn_posterior.keys())

    # Print summary statistics and get per-chromosome aneuploidy calls
    aneuploid_chromosomes = summary_statistics(data, map_estimates, cn_posterior, prob_threshold=args.prob_threshold)

    # Separate samples into normal and aneuploid groups
    normal_samples = []
    aneuploid_samples_list = []

    for sample_idx in range(data.n_samples):
        if len(aneuploid_chromosomes[sample_idx]) == 0:
            normal_samples.append(sample_idx)
        else:
            aneuploid_samples_list.append(sample_idx)

    print(f"\n{'='*80}")
    print("GENERATING COMBINED RESULTS PLOTS")
    print(f"{'='*80}")
    print(f"Normal samples: {len(normal_samples)}")
    print(f"Aneuploid samples: {len(aneuploid_samples_list)}")

    # Generate bin variance/bias plot once (global across all samples)
    print("\nGenerating global bin variance and bias plot...")
    plot_bin_variance_bias(data, map_estimates, output_dir=output_dir)

    # Generate combined plots for aneuploid samples
    n_aneuploid_to_plot = len(aneuploid_samples_list)
    print(f"\nGenerating plots for {n_aneuploid_to_plot} aneuploid samples...")
    for i, sample_idx in enumerate(aneuploid_samples_list[:n_aneuploid_to_plot]):
        sample_name = data.sample_ids[sample_idx]
        aneuploid_chrs = aneuploid_chromosomes[sample_idx]
        chr_list = ', '.join([f"{chr_name}(CN={cn})" for chr_name, cn, _ in aneuploid_chrs])
        print(f"  {i+1}/{n_aneuploid_to_plot}: {sample_name} (sample_idx={sample_idx}) - {chr_list}")
        plot_sample(data, map_estimates, cn_posterior, output_dir=output_dir, sample_idx=sample_idx, is_aneuploid=True, aneuploid_chrs=aneuploid_chrs)

    # Generate combined plots for normal samples
    n_normal_to_plot = len(normal_samples)
    print(f"\nGenerating plots for {n_normal_to_plot} normal samples...")
    for i, sample_idx in enumerate(normal_samples[:n_normal_to_plot]):
        sample_name = data.sample_ids[sample_idx]
        print(f"  {i+1}/{n_normal_to_plot}: {sample_name} (sample_idx={sample_idx})")
        plot_sample(data, map_estimates, cn_posterior, output_dir=output_dir, sample_idx=sample_idx, is_aneuploid=False, aneuploid_chrs=None)

    print(f"{'='*80}\n")

    # Create results dataframe with copy number calls and probabilities
    results = []
    cn_probs = cn_posterior['cn_posterior']

    for i in range(data.n_bins):
        for j in range(data.n_samples):
            cn_map = map_estimates['cn'][i, j]
            cn_prob = cn_probs[i, j, :]
            
            result = {
                'chr': data.chr[i],
                'start': data.start[i],
                'end': data.end[i],
                'sample': data.sample_ids[j],
                'observed_depth': data.depth[i, j].cpu().numpy(),
                'cn_map': int(cn_map),
                'cn_prob_0': cn_prob[0],
                'cn_prob_1': cn_prob[1],
                'cn_prob_2': cn_prob[2],
                'cn_prob_3': cn_prob[3],
                'max_prob': cn_prob.max(),
                'bin_bias': map_estimates['bin_bias'][i],
                'bin_var': map_estimates['bin_var'][i],
                'sample_var': map_estimates['sample_var'][j],
            }
            results.append(result)

    results_df = pd.DataFrame(results)

    # Save results to TSV (gzipped)
    output_file = os.path.join(output_dir, 'bin_stats.tsv.gz')
    results_df.to_csv(output_file, sep='\t', index=False, compression='gzip')
    print(f"Results saved to {output_file}")

    # Save per-chromosome statistics for all samples (both normal and aneuploid)
    chromosome_stats = []
    cn_probs = cn_posterior['cn_posterior']
    
    # Get unique chromosomes
    unique_chrs = np.unique(data.chr)
    
    for sample_idx in range(data.n_samples):
        sample_name = data.sample_ids[sample_idx]
        
        for chr_name in unique_chrs:
            # Get chromosome mask and bins
            chr_mask = data.chr == chr_name
            chr_bins_count = np.sum(chr_mask)
            
            if chr_bins_count == 0:
                continue
            
            # Get CN calls for this chromosome
            chr_cn = map_estimates['cn'][chr_mask, sample_idx]
            cn_counts = np.bincount(chr_cn, minlength=6)
            most_common_cn = np.argmax(cn_counts)
            
            # Calculate mean probability of the called CN state
            chr_cn_probs = cn_probs[chr_mask, sample_idx, most_common_cn]
            mean_prob = chr_cn_probs.mean()
            
            # Get depth statistics
            chr_depths = data.depth[chr_mask, sample_idx].cpu().numpy()
            
            # Determine if this is an aneuploidy
            is_aneuploid = False
            is_high_confidence = mean_prob > 0.9  # Using same threshold as summary_statistics
            
            # Check if it's in the aneuploid_chromosomes dictionary
            if sample_idx in aneuploid_chromosomes:
                for aneu_chr, aneu_cn, aneu_prob in aneuploid_chromosomes[sample_idx]:
                    if aneu_chr == chr_name:
                        is_aneuploid = True
                        break
            
            chromosome_stats.append({
                'sample': sample_name,
                'chromosome': chr_name,
                'copy_number': int(most_common_cn),
                'mean_cn_probability': float(mean_prob),
                'is_aneuploid': is_aneuploid,
                'is_high_confidence': is_high_confidence,
                'n_bins': int(chr_bins_count),
                'mean_depth': float(np.mean(chr_depths)),
                'median_depth': float(np.median(chr_depths)),
                'std_depth': float(np.std(chr_depths)),
                'sample_var': float(map_estimates['sample_var'][sample_idx]),
            })
    
    # Create and save chromosome stats dataframe
    chromosome_stats_df = pd.DataFrame(chromosome_stats)
    chromosome_stats_df = chromosome_stats_df.sort_values(['sample', 'chromosome'])
    
    chromosome_stats_output = os.path.join(output_dir, 'chromosome_stats.tsv')
    chromosome_stats_df.to_csv(chromosome_stats_output, sep='\t', index=False)
    print(f"\nChromosome statistics saved to {chromosome_stats_output}")


if __name__ == "__main__":
    main()
