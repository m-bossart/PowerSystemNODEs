import tensorflow as tf
import sys
from sklearn.naive_bayes import GaussianNB, MultinomialNB, ComplementNB, BernoulliNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn import svm

def NaiiveBayes(x_train, y_train, NB_type):
    if NB_type == 'Gaussian':
        model = GaussianNB()
    elif NB_type == 'Multinomial':
        model = MultinomialNB()
    elif NB_type == 'Complement':
        model = ComplementNB()
    elif NB_type == 'Bernoulli':
        model = BernoulliNB()
    else:
        raise ValueError('Undefined {} type for NB'.format(NB_type))
    model.fit(x_train, y_train)
    return model


def KNN(x_train, y_train, k):
    model = KNeighborsClassifier(n_neighbors=k)
    return model


def LR(x_train, y_train, c):
    model = LogisticRegression(random_state=0).fit(x_train, y_train)
    return model

def SVM(x_train, y_train, c, kernel):
    model = svm.SVC(C=c, kernel=kernel)
    model.fit(x_train, y_train)
    return model

def lstm(hparams, input_size, network_config):
    inputs = tf.keras.Input(shape=(input_size, 1))

    prev = inputs
    layers = network_config["layers"]
    units = [i[0] for i in layers]
    directions = [i[1] for i in layers]

    for i, (size, direction) in enumerate(zip(units, directions)):
        return_sequences = False if i >= len(units) - 1 else True
        if direction == 1:
            prev = tf.keras.layers.LSTM(size, return_sequences=return_sequences)(prev)
        elif direction == 2:
            prev = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(
                size, return_sequences=return_sequences))(prev)
        else:
            sys.exit("{} is not a valid LSTM direction".format(direction))
    outputs = tf.keras.layers.Dense(hparams['num_classes'],
                                    activation=hparams['output_activation'])(prev)

    return inputs, outputs


def cnn(hparams, input_size, network_config):
    inputs = tf.keras.layers.Input(tuple(network_config['input_shape']))
    layers = network_config['layers']
    prev = inputs

    for i in layers:
        prev = tf.keras.layers.Conv2D(i[0], i[1],
                                      activation=network_config["activation"])(prev)
    shared = tf.keras.layers.GlobalMaxPooling2D()(prev)
    shared = tf.keras.layers.Flatten()(prev)
    outputs = tf.keras.layers.Dense(hparams['num_classes'],
                                    activation=hparams['output_activation'])(shared)
    return inputs, outputs


def fully_connected(hparams, input_size, network_config):
    inputs = tf.keras.Input(shape=(input_size,), dtype=tf.dtypes.float32)
    layers = build_dense_layers(inputs, network_config, name="layers")
    outputs = tf.keras.layers.Dense(hparams['num_classes'],
                                    activation=hparams['output_activation'])(layers)
    return inputs, outputs


def build_dense_layers(inputs, network_config, name):
    hidden_nodes = network_config["hidden_nodes"]
    prev = inputs

    for i, size in enumerate(hidden_nodes):
        layer = tf.keras.layers.Dense(
            size, activation=network_config['activation'], name="{}_{}".format(name, i))(prev)
        prev = tf.keras.layers.Dropout(network_config['dropout'])(layer)

    return prev
