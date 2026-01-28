#!/bin/python

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


def plot_combined_results(data: AneuploidyData, map_estimates: dict, cn_posterior: dict, output_dir: str,
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
    fig, axes = plt.subplots(3, 1, figsize=(16, 12))
    
    # Extract data
    observed = data.depth[:, sample_idx].cpu().numpy()
    cn = map_estimates['cn'][:, sample_idx]
    sample_var_all = map_estimates['sample_var'].flatten()  # All samples' variance factors
    sample_var_current = sample_var_all[sample_idx]  # Current sample's variance factor
    cn_probs = cn_posterior['cn_posterior'][:, sample_idx, :]  # (n_bins, n_states)
    
    colors = [ '#004D40', '#FFC107',  '#1E88E5', '#D81B60', '#38006B']
    x = np.arange(len(observed))
    
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
    ax.plot(x_transformed, cn, 'o', alpha=0.5, markersize=6,
            label=f'MAP copy state', color='black')
    # Then plot observed depth on top so it's visible
    ax.plot(x_transformed, observed, 'r-', linewidth=0.8, alpha=0.7, label='Normalized read depth')
    ax.set_title(f'Sample: {data.sample_ids[sample_idx]}')
    ax.legend(loc='lower left')
    ax.grid(True, axis='y', alpha=1, linestyle='-', linewidth=1)
    ax.set_ylim([-0.5, 4.5])
    ax.set_xlim([x_transformed.min(), x_transformed.max()])
    
    # Plot 2: Stacked area plot of CN probabilities
    ax = axes[1]
    ax.stackplot(x_transformed, cn_probs.T, labels=[f'CN={i}' for i in range(6)], 
                alpha=0.7, colors=colors)
    ax.set_ylabel('Copy Number Probability')
    ax.legend(loc='upper left', ncol=6)
    ax.set_ylim([0, 1])
    ax.set_xlim([x_transformed.min(), x_transformed.max()])
    
    # Plot 3: Sample variance distribution
    ax = axes[2]
    ax.hist(np.sqrt(sample_var_all), bins=30, alpha=0.7, edgecolor='black', linewidth=0.5, color='gray')
    ax.axvline(np.sqrt(sample_var_current), color='red', linestyle='--', linewidth=2, 
              label=f'This sample: {np.sqrt(sample_var_current):.3f}')
    ax.set_xlabel('MAP sample std')
    ax.set_ylabel('Count')
    ax.set_yscale('log')
    ax.legend(loc='upper right')
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
    aneu_suffix = '_ANEU' if is_aneuploid else ''
    filename = f'combined_results_{sample_name}{aneu_suffix}.png'
    output_path = os.path.join(output_dir, filename)
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
    ax.plot(x, np.sqrt(bin_var), 'o', alpha=0.5, markersize=3, color='purple', label='Bin variance')
    ax.set_ylabel('Bin std')
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


def read_data(pattern: str) -> pd.DataFrame:
    files = glob.glob(pattern, recursive=True)
    
    print(f"Found {len(files)} files matching pattern: {pattern}")
    
    # Load first file to establish the base dataframe
    df_merged = None
    
    # Track successful and failed files
    successful_files = []
    failed_files = []
    
    # Load each file
    for file_path in files:
        print(f"Loading: {file_path}")
        try:
            # Read gzipped TSV file
            df_file = pd.read_csv(file_path, sep='\t', compression='gzip')
            
            # Get source file identifier
            source_file = str(Path(file_path).parent.parent.name)
            
            # Rename chromosome column if it exists
            if '#Chr' in df_file.columns:
                df_file['Chr'] = df_file['#Chr']
                df_file = df_file.drop('#Chr', axis=1)
            
            # Create bin identifier
            df_file['Bin'] = df_file['Chr'].astype(str) + ':' + df_file['Start'].astype(str) + '-' + df_file['End'].astype(str)
            df_file = df_file.set_index('Bin')
            
            # Get sample columns (exclude metadata)
            metadata_cols = ['Chr', 'Start', 'End']
            sample_cols = [col for col in df_file.columns if col not in metadata_cols]
            
            if df_merged is None:
                # First file: keep metadata + samples
                df_merged = df_file.copy()
                df_merged['source_file'] = source_file
                print(f"  Initialized with {len(df_merged)} bins and {len(sample_cols)} samples")
            else:
                # Subsequent files: merge sample columns only
                # Keep only sample columns from this file
                df_samples = df_file[sample_cols].copy()
                
                # Check for duplicate sample columns and drop them
                existing_sample_cols = [col for col in df_merged.columns if col not in ['Chr', 'Start', 'End', 'source_file']]
                duplicate_samples = [col for col in sample_cols if col in existing_sample_cols]
                
                if duplicate_samples:
                    print(f"  WARNING: Found {len(duplicate_samples)} duplicate sample(s) in {file_path}")
                    print(f"  Dropping duplicate samples: {duplicate_samples[:5]}{'...' if len(duplicate_samples) > 5 else ''}")
                    df_samples = df_samples.drop(columns=duplicate_samples)
                    sample_cols = [col for col in sample_cols if col not in duplicate_samples]
                
                # Skip this file if no new samples remain after dropping duplicates
                if len(sample_cols) == 0:
                    print(f"  Skipping file - all samples already exist in merged dataframe")
                    continue
                
                # Validate that bins match before merging
                existing_bins = set(df_merged.index)
                new_bins = set(df_samples.index)
                
                if existing_bins != new_bins:
                    # Bins don't match - try coordinate adjustments before failing
                    print(f"\n  WARNING: Bin mismatch detected for {file_path}")
                    print(f"  Attempting coordinate adjustment...")
                    
                    # Try adding 1 to coordinates
                    df_file_adjusted = df_file.copy()
                    df_file_adjusted['Start'] = df_file_adjusted['Start'] + 1
                    df_file_adjusted['End'] = df_file_adjusted['End'] + 1
                    df_file_adjusted['Bin'] = (df_file_adjusted['Chr'].astype(str) + ':' + 
                                               df_file_adjusted['Start'].astype(str) + '-' + 
                                               df_file_adjusted['End'].astype(str))
                    df_file_adjusted = df_file_adjusted.set_index('Bin')
                    df_samples_adjusted = df_file_adjusted[sample_cols].copy()
                    new_bins_adjusted = set(df_samples_adjusted.index)
                    
                    if existing_bins == new_bins_adjusted:
                        print(f"  ✓ Success: Coordinates adjusted by +1")
                        df_samples = df_samples_adjusted
                    else:
                        # Try subtracting 1 instead
                        df_file_adjusted = df_file.copy()
                        df_file_adjusted['Start'] = df_file_adjusted['Start'] - 1
                        df_file_adjusted['End'] = df_file_adjusted['End'] - 1
                        df_file_adjusted['Bin'] = (df_file_adjusted['Chr'].astype(str) + ':' + 
                                                   df_file_adjusted['Start'].astype(str) + '-' + 
                                                   df_file_adjusted['End'].astype(str))
                        df_file_adjusted = df_file_adjusted.set_index('Bin')
                        df_samples_adjusted = df_file_adjusted[sample_cols].copy()
                        new_bins_adjusted = set(df_samples_adjusted.index)
                        
                        if existing_bins == new_bins_adjusted:
                            print(f"  ✓ Success: Coordinates adjusted by -1")
                            df_samples = df_samples_adjusted
                        else:
                            # Neither adjustment worked - show error and raise
                            only_in_existing = existing_bins - new_bins
                            only_in_new = new_bins - existing_bins
                            
                            print("\n" + "="*80)
                            print("ERROR: Bin mismatch detected!")
                            print("="*80)
                            print(f"\nFile: {file_path}")
                            print(f"Existing dataframe has {len(existing_bins)} bins")
                            print(f"New file has {len(new_bins)} bins")
                            print(f"Bins only in existing: {len(only_in_existing)}")
                            print(f"Bins only in new file: {len(only_in_new)}")
                            
                            if only_in_existing:
                                print(f"\nExample bins only in existing dataframe (first 5):")
                                for bin_id in list(only_in_existing)[:5]:
                                    print(f"  {bin_id}")
                            
                            if only_in_new:
                                print(f"\nExample bins only in new file (first 5):")
                                for bin_id in list(only_in_new)[:5]:
                                    print(f"  {bin_id}")
                            
                            raise ValueError(
                                f"Bin mismatch: Cannot merge {file_path} because bins don't align. "
                                f"Expected {len(existing_bins)} bins, found {len(new_bins)} bins. "
                                f"Coordinate adjustments (±1) did not resolve the mismatch."
                            )
                
                # Merge with existing dataframe
                df_merged = df_merged.join(df_samples, how='outer')
                print(f"  Merged {len(sample_cols)} samples, total bins: {len(df_merged)}")
            
            # Mark as successful if we got here
            successful_files.append(file_path)
                
        except Exception as e:
            print(f"Error loading {file_path}: {e}")
            failed_files.append((file_path, str(e)))
    
    if df_merged is not None:
        # Remove duplicate bins (keep first occurrence)
        df_merged = df_merged[~df_merged.index.duplicated(keep='first')]
        
        print(f"\nCombined dataframe shape: {df_merged.shape}")
        print(f"Total bins: {len(df_merged)}")
        print(f"Total columns: {len(df_merged.columns)}")
        
        # Count sample columns
        metadata_cols = ['Chr', 'Start', 'End', 'source_file']
        sample_count = len([col for col in df_merged.columns if col not in metadata_cols])
        print(f"Sample columns: {sample_count}")
    else:
        print("No files were successfully loaded")
        df_merged = pd.DataFrame()
    
    # Print consolidated file loading summary
    print("\n" + "="*80)
    print("FILE LOADING SUMMARY")
    print("="*80)
    print(f"Total files found: {len(files)}")
    print(f"Successfully loaded: {len(successful_files)}")
    print(f"Failed to load: {len(failed_files)}")
    
    if failed_files:
        print("\nFailed files:")
        for file_path, error in failed_files:
            print(f"  - {file_path}")
            print(f"    Error: {error}")
    print("="*80 + "\n")
    
    return df_merged


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


def main():
    """
    Main function to run aneuploidy detection pipeline.
    """
    # Find all matching files
    bins_pattern = "/Users/markw/Work/talkowski/sv-pipe-testing/mw_ploidy/cmg/data/m2_batch32_ploidy_plots/ploidy_est/binwise_estimated_copy_numbers.bed.gz"
    output_dir = "/Users/markw/Work/talkowski/sv-pipe-testing/mw_ploidy/cmg/aneuploidy_pyro_output"

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Delete all existing combined_results_*.png files
    combined_results_pattern = os.path.join(output_dir, "combined_results_*.png")
    existing_files = glob.glob(combined_results_pattern)
    if existing_files:
        print(f"Deleting {len(existing_files)} existing combined_results_*.png files...")
        for file_path in existing_files:
            os.remove(file_path)
        print("Cleanup complete.")

    df = read_data(bins_pattern)
    # TODO filter bins
    #df = df[df['Chr'] == 'chrY']
    #df = df[df['Chr'].isin({'chrX', 'chrY'})]
    #df = df[df['Chr'].isin({'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'})]
    #df = df.loc[:, ['Chr', 'Start', 'End', 'source_file', '__cdh1355__1a8a28'] + [x for x in df.columns[2:12].values.tolist() if x != '__cdh1355__1a8a28']]
    #df = df.loc[:, ['Chr', 'Start', 'End', 'source_file'] + df.columns[2:152].values.tolist()]

    # Filter low quality bins
    df = filter_low_quality_bins(df,
                                autosome_median_min=1.75,
                                autosome_median_max=2.25,
                                autosome_mad_max=0.25,
                                chrX_median_min=0.75,
                                chrX_median_max=2.25,
                                chrX_mad_max=0.5,
                                chrY_median_min=0.0,
                                chrY_median_max=1.25,
                                chrY_mad_max=0.5)

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
    device = 'cpu'  # Use 'cuda' if GPU available
    dtype = torch.float32
    inference_method = 'svi'

    data = AneuploidyData(df_subset, device=device, dtype=dtype)

    # Initialize model
    model = AneuploidyModel(
        n_states=4,           # CN states: 0, 1, 2, 3
        alpha_ref=1.0,       # Strong prior on CN=2; 100 for autosome, 1 for allosome
        alpha_non_ref=1.0,    # Weak prior on other states
        var_bias_bin=0.01,     # Moderate bin-to-bin bias variation
        var_sample=0.01,       # Moderate sample variance
        var_bin=0.01,          # Moderate bin variance
        device=device,
        dtype=dtype,
        guide_type='delta'      # MAP inference
    )

    # Train model
    if inference_method == 'svi':
        model.train(
            data=data,
            max_iter=1000,         # Number of training iterations (increase for better convergence)
            lr_init=0.01,         # Initial learning rate
            lr_min=0.01,         # Minimum learning rate
            lr_decay=10000,         # Learning rate decay constant
            log_freq=50,          # Log every 50 iterations
            jit=True             # Set to True for faster training (may have compatibility issues)
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
    cn_posterior = model.run_discrete_inference(data, n_samples=1000, log_freq=100)
    print("cn_posterior keys:", cn_posterior.keys())

    # Print summary statistics and get per-chromosome aneuploidy calls
    aneuploid_chromosomes = summary_statistics(data, map_estimates, cn_posterior)

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

    # Generate combined plots for up to 10 normal samples
    n_normal_to_plot = min(10, len(normal_samples))
    print(f"\nGenerating plots for {n_normal_to_plot} normal samples...")
    for i, sample_idx in enumerate(normal_samples[:n_normal_to_plot]):
        sample_name = data.sample_ids[sample_idx]
        print(f"  {i+1}/{n_normal_to_plot}: {sample_name} (sample_idx={sample_idx})")
        plot_combined_results(data, map_estimates, cn_posterior, output_dir=output_dir, sample_idx=sample_idx, is_aneuploid=False, aneuploid_chrs=None)

    # Generate combined plots for up to 10 aneuploid samples
    n_aneuploid_to_plot = min(10, len(aneuploid_samples_list))
    print(f"\nGenerating plots for {n_aneuploid_to_plot} aneuploid samples...")
    for i, sample_idx in enumerate(aneuploid_samples_list[:n_aneuploid_to_plot]):
        sample_name = data.sample_ids[sample_idx]
        aneuploid_chrs = aneuploid_chromosomes[sample_idx]
        chr_list = ', '.join([f"{chr_name}(CN={cn})" for chr_name, cn, _ in aneuploid_chrs])
        print(f"  {i+1}/{n_aneuploid_to_plot}: {sample_name} (sample_idx={sample_idx}) - {chr_list}")
        plot_combined_results(data, map_estimates, cn_posterior, output_dir=output_dir, sample_idx=sample_idx, is_aneuploid=True, aneuploid_chrs=aneuploid_chrs)

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
            }
            results.append(result)

    results_df = pd.DataFrame(results)
    print(f"Results dataframe shape: {results_df.shape}")
    print(results_df.head(20))

    # Save results to TSV
    output_file = os.path.join(output_dir, 'aneuploidy_results.tsv')
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to {output_file}")

    # Save per-chromosome high-confidence aneuploidy calls
    chromosome_aneuploidies = []
    for sample_idx, aneuploid_chrs in aneuploid_chromosomes.items():
        sample_name = data.sample_ids[sample_idx]
        for chr_name, cn, mean_prob in aneuploid_chrs:
            # Get chromosome bins for additional info
            chr_mask = data.chr == chr_name
            chr_bins_data = results_df[(results_df['chr'] == chr_name) & (results_df['sample'] == sample_name)]
            
            # Convert observed_depth arrays to scalars for statistics
            depths = np.array([float(d) if hasattr(d, '__iter__') and not isinstance(d, str) else float(d) 
                            for d in chr_bins_data['observed_depth']])
            
            chromosome_aneuploidies.append({
                'sample': sample_name,
                'chromosome': chr_name,
                'copy_number': int(cn),
                'mean_cn_probability': mean_prob,
                'n_bins': chr_bins_data.shape[0],
                'mean_depth': float(np.mean(depths)),
                'median_depth': float(np.median(depths)),
                'std_depth': float(np.std(depths)),
            })

    if chromosome_aneuploidies:
        aneuploidy_df = pd.DataFrame(chromosome_aneuploidies)
        aneuploidy_df = aneuploidy_df.sort_values(['sample', 'chromosome'])
        print(f"Total chromosomal aneuploidies: {len(aneuploidy_df)}")
    else:
        print("No high-confidence chromosomal aneuploidies detected")
        # Create empty dataframe with expected columns
        aneuploidy_df = pd.DataFrame(columns=['sample', 'chromosome', 'copy_number', 'mean_cn_probability', 
                                            'n_bins', 'mean_depth', 'median_depth', 'std_depth'])

    aneuploidy_output = os.path.join(output_dir, 'high_confidence_aneuploidies.tsv')
    aneuploidy_df.to_csv(aneuploidy_output, sep='\t', index=False)
    print(f"High-confidence aneuploidy calls saved to {aneuploidy_output}")


if __name__ == "__main__":
    main()
