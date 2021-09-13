import os
from collections import defaultdict

from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import numpy as np
import tensorflow as tf

def evaluate_model(model, x_test, y_test, train_type):
    y_pred = model.predict(x_test)

    print('Mean squared error: %.2f' % mean_squared_error(y_test, y_pred, squared=False))
    print('Integer - Mean squared error: %.2f \n' % mean_squared_error(y_test, list(map(int, y_pred)), squared=False))
#    import pdb;pdb.set_trace()
    return mean_squared_error(y_test, y_pred)

def preprocess_data(x_train, x_test, preprocess):

    if preprocess == "scalar":
        scaler = StandardScaler()
    elif preprocess == "minmax":
        scaler = MinMaxScaler()
    else:
        raise ValueError("preprocessing method not available")

    train_preprocessed = scaler.fit_transform(x_train)

    test_preprocessed = scaler.transform(x_test)

    return train_preprocessed, test_preprocessed
    


def build_parameter_combinations(params):
    base = [{}]
    for param, values in params.items():
        frontier = []
        for b in base:
            for value in values:
                c = {k: v for k, v in b.items()}
                c[param] = value
                frontier.append(c)
        base = frontier

    return frontier
