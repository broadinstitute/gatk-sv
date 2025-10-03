import argparse
import json
import logging
import numpy as np
import os
import pyro
from pyro.infer import JitTraceEnum_ELBO, TraceEnum_ELBO
from pyro.infer.svi import SVI
import torch

from .cnv_model import SVDepthData, SVDepthPyroModel
from . import cnv_io


def train_epoch(svi: SVI, data: SVDepthData, epoch: int, scheduler: pyro.optim.LambdaLR) -> float:
    loss = svi.step(counts=data.counts, bin_size=data.bin_size, sample_ploidy=data.sample_ploidy,
                    sample_depth=data.sample_per_base_depth)
    scheduler.step(epoch)
    return loss


def create_scheduler(lr_min: float, lr_init: float, lr_decay: float, beta1: float, beta2: float):
    return pyro.optim.LambdaLR({
        'optimizer': torch.optim.Adam,
        'optim_args': {'lr': 1., 'betas': (beta1, beta2)},
        'lr_lambda': lambda k: lr_min + (lr_init - lr_min) * np.exp(-k / lr_decay)
    })


def run_training(model: SVDepthPyroModel,
                 data: SVDepthData,
                 max_iter: int,
                 lr_min: float,
                 lr_init: float,
                 lr_decay: float,
                 adam_beta1: float,
                 adam_beta2: float,
                 iter_log_freq: int,
                 jit: bool):
    logging.info("Initializing model...")
    pyro.clear_param_store()
    scheduler = create_scheduler(lr_min=lr_min, lr_init=lr_init, lr_decay=lr_decay,
                                 beta1=adam_beta1, beta2=adam_beta2)
    if jit:
        loss = JitTraceEnum_ELBO()
    else:
        loss = TraceEnum_ELBO()
    svi = SVI(model.model, model.guide, optim=scheduler, loss=loss)
    logging.info("Running model training...")
    # Run training loop.  Use try to allow for keyboard interrupt.
    try:
        for epoch in range(max_iter):
            # Train, and keep track of training loss.
            total_epoch_loss = train_epoch(svi=svi, data=data, epoch=epoch, scheduler=scheduler)
            model.loss['epoch'].append(epoch)
            model.loss['elbo'].append(-total_epoch_loss)
            if (epoch + 1) % iter_log_freq == 0:
                logging.info("[epoch %04d]  training loss: %.4f" % (epoch + 1, total_epoch_loss))
        logging.info("Training procedure complete after {:d} epochs, loss: {:.4f}".format(max_iter, total_epoch_loss))

    # Exception allows program to continue after ending training prematurely.
    except KeyboardInterrupt:
        logging.info("Training procedure stopped by keyboard interrupt.")


def run(args: dict, default_dtype: torch.dtype = torch.float32):
    base_path = os.path.join(args['output_dir'], args['output_name'])
    log_path = base_path + ".log.txt"
    logging.basicConfig(filename=log_path,
                        filemode='w',
                        level=logging.INFO,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    pyro.enable_validation(True)
    pyro.distributions.enable_validation(True)
    pyro.clear_param_store()

    np.random.seed(args['random_seed'])
    torch.random.manual_seed(args['random_seed'])
    pyro.set_rng_seed(args['random_seed'])

    model = SVDepthPyroModel(k=args['num_states'], tensor_dtype=default_dtype,  read_length=args['read_length'],
                             var_phi_bin=args['var_phi_bin'], var_phi_sample=args['var_phi_sample'],
                             alpha_ref=args['alpha_ref'], alpha_non_ref=args['alpha_non_ref'],
                             mu_eps=args['mu_eps'], device=args['device'])
    data = cnv_io.load_data(variants_file_path=args['data_file'], mean_coverage_path=args['sample_depth_file'],
                            samples_path=args['samples_file'], tensor_dtype=default_dtype, device=args['device'])
    if data is None:
        logging.warning("No intervals found.")
    else:
        logging.info("Training with {:d} intervals and {:d} samples...".format(data.counts.shape[0], data.counts.shape[1]))
        run_training(model=model, data=data, max_iter=args['max_iter'], lr_min=args['lr_min'], lr_init=args['lr_init'],
                     lr_decay=args['lr_decay'], adam_beta1=args['adam_beta1'], adam_beta2=args['adam_beta2'],
                     iter_log_freq=args['iter_log_freq'], jit=args['jit'])
        save_data(base_path=base_path, model=model, data=data, contigs=data.contigs, sample_ids=data.samples)


def save_data(base_path: str, model: SVDepthPyroModel, data: SVDepthData, contigs: list, sample_ids: list):
    cnv_io.save_tensors(data=data, base_path=base_path)
    cnv_io.save_list(data=contigs, path=base_path + ".contigs.list")
    cnv_io.save_list(data=sample_ids, path=base_path + ".sample_ids.list")
    save_model(model=model, base_path=base_path)


def save_model(model: SVDepthPyroModel, base_path: str):
    with open(base_path + '.vars.json', 'w') as f:
        json.dump({
            'k': model.k,
            'var_phi_bin': model.var_phi_bin,
            'var_phi_sample': model.var_phi_sample,
            'mu_eps': model.mu_eps,
            'alpha_ref': model.alpha_ref,
            'alpha_non_ref': model.alpha_non_ref,
            'read_length': model.read_length,
            'loss': model.loss,
            'guide_loc': pyro.param('AutoDiagonalNormal.loc').tolist(),
            'guide_scale': pyro.param('AutoDiagonalNormal.scale').tolist()
        }, f)
    pyro.get_param_store().save(base_path + '.param_store.pyro')


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svgenotyper cnv-train',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--data_file", type=str, required=True)
    parser.add_argument("--sample_depth_file", type=str, required=True)
    parser.add_argument("--samples_file", type=str, required=True)
    parser.add_argument("--output_name", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)

    parser.add_argument("--device", type=str, default='cpu')
    parser.add_argument("--jit", action="store_true")

    parser.add_argument("--num_states", type=int, default=5)
    parser.add_argument("--random_seed", type=int, default=92837488)
    parser.add_argument("--read_length", type=int, default=150)
    parser.add_argument("--max_iter", type=int, default=2000)
    parser.add_argument("--iter_log_freq", type=int, default=50)

    parser.add_argument("--var_phi_bin", type=float, default=0.1)
    parser.add_argument("--var_phi_sample", type=float, default=0.1)
    parser.add_argument("--alpha_ref", type=float, default=18)
    parser.add_argument("--alpha_non_ref", type=float, default=1)
    parser.add_argument("--mu_eps", type=float, default=.01)

    parser.add_argument("--lr_decay", type=float, default=1000)
    parser.add_argument("--lr_min", type=float, default=0.001)
    parser.add_argument("--lr_init", type=float, default=0.01)
    parser.add_argument("--adam_beta1", type=float, default=0.9)
    parser.add_argument("--adam_beta2", type=float, default=0.999)

    args = parser.parse_args(argv)

    print('Arguments are', args)
    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    run(args)


if __name__ == '__main__':
    main()