import argparse
import traceback
import yaml

from src.utils import build_parameter_combinations
from src.model import train

import os
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"]="1,2,3"

def hyper_tuning(hparams, tuning_path):
    with open(tuning_path) as f:
        tuning_params = yaml.load(f, Loader=yaml.FullLoader)

    combinations = build_parameter_combinations(tuning_params)
    for combination_id, combination in enumerate(combinations):
        print("Combination: {} / {}".format(combination_id, len(combinations)))
        for key, value in combination.items():
            node = hparams
            path = key.split("|")
            for i, k in enumerate(path):
                if i == len(path) - 1:
                    node[k] = value
                else:
                    node = node[k]

        try:
            train(hparams)
        except:
            traceback.print_exc()


if __name__ == "__main__":
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("-c", "--config", type=str, default='config.yml',
                             help="Experiment Configurations ")
    args_parser.add_argument("-t", "--tuning", type=str, default=None,
                             help="Experiment Configurations ")
    run_args = args_parser.parse_args()
    config_path = run_args.config
    tuning_path = run_args.tuning

    with open(config_path) as f:
        hparams = yaml.load(f, Loader=yaml.FullLoader)

    if tuning_path:
        hyper_tuning(hparams, tuning_path)
    else:
        train(hparams)
