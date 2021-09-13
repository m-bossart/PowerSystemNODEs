from functools import reduce

import tensorflow as tf
from tensorflow.keras.callbacks import Callback
from tensorflow.keras.models import clone_model

from utils import evaluate_model


class EvaluateCallback(Callback):
    def __init__(self, data_generator, name, eval_period, hparams):
        self.data_generator = data_generator
        self.name = name
        self.eval_period = eval_period
        self.best_models = {"sum": {"epoch": 0, "score": 0, "results": None, "model": None},
                            "geo": {"epoch": 0, "score": 0, "results": None, "model": None}}

        self.current_epoch = 0
        self.hparams = hparams

    def on_epoch_end(self, epoch, logs=None):
        self.current_epoch += 1
        if (epoch + 1) % self.eval_period == 0:
            results = evaluate_model(self.model, self.data_generator, self.hparams)
            f1_eval = {}

            f1_eval["sum"] = sum([results[x+"_f1"] for x in self.data_generator.selected_tasks])
            f1_eval["geo"] = reduce(lambda x, y: x * y,
                                    [results[x+"_f1"] for x in self.data_generator.selected_tasks])

            for k in f1_eval:
                if f1_eval[k] > self.best_models[k]["score"]:
                    print("\n", "*** New Best {} {} Evaluation: ".format(k, self.name), results)
                    self.best_models[k]["results"] = results
                    self.best_models[k]["score"] = f1_eval[k]
                    self.best_models[k]["model"] = clone_model(self.model)
                    self.best_models[k]["model"].set_weights(self.model.get_weights())
                    self.best_models[k]["epoch"] = self.current_epoch


class LearningRateCallback():
    def __init__(self, hparams):
        self.decay_after = hparams['decay_after']
        self.base = hparams['base']
        self.decay = hparams['decay']
        self.total_epochs = 0

    def calculate(self, epoch):
        self.total_epochs += 1
        if self.total_epochs < self.decay_after:
            return self.base
        else:
            return self.base * tf.math.exp(
                self.decay * (self.decay_after - self.total_epochs))
