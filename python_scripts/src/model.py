import pickle
import tensorflow as tf

from tensorflow.keras.callbacks import LearningRateScheduler
from keras import backend as K
from networks import lstm, cnn, fully_connected, NaiiveBayes, KNN, LR, SVM
from utils import evaluate_model, preprocess_data
from data_utils import get_data

def train(hparams):
    input_size = hparams["input_size"]

    (train_data, Y_dev), (test_data, Y_test), (all_data, Y), features_str = get_data(hparams["data_path"])
    train_data, test_data = preprocess_data(train_data, test_data, hparams["preprocess"])

    if hparams['training_type'] == 'DL':

        if hparams["network_type"] == "lstm":
            inputs, outputs = lstm(hparams,
                                   input_size,
                                   hparams["network_config"]["lstm"])
        elif hparams["network_type"] == "cnn":
            inputs, outputs = cnn(hparams,
                                  input_size,
                                  hparams["network_config"]["cnn"])
        elif hparams["network_type"] == "fully_connected":
            inputs, outputs = fully_connected(hparams,
                                              input_size,
                                              hparams["network_config"]["fully_connected"])
        else:
            raise ValueError('Undefined {} network for DL'.format(hparams["network_type"]))

        model = tf.keras.Model(inputs=inputs, outputs=outputs, name="RUL_estimator")
        model.summary()

        model.compile(optimizer=tf.keras.optimizers.Adam(hparams["learning_rate"]),
                      loss=hparams["loss"])
        model.fit(train_data, Y_dev, verbose=1, shuffle=True, validation_split=0.1,
                    epochs=hparams["epochs"], batch_size=hparams["batch_size"])


    elif hparams["training_type"] == 'ML':

        if hparams["network_type"] == 'NB':
            model = NaiiveBayes(x_train, y_train, hparams["network_config"]["NB"]["type"])
        elif hparams["network_type"] == 'KNN':
            model = KNN(x_train, y_train, hparams["network_config"]["KNN"]["k"])
        elif hparams["network_type"] == 'logistic_regression':
            model = LR(x_train, y_train, hparams["network_config"]["logistic_regression"]["c"])
        elif hparams["network_type"] == 'SVM':
            model = SVM(x_train, y_train, hparams["network_config"]["SVM"]["c"],
                                          hparams["network_config"]["SVM"]["kernel"])
        else:
            raise ValueError('Undefined {} network type for ML'.format(hparams["network_type"]))

    if hparams["network_type"] == "cnn":
        model.evaluate(ds_test)
    else:
        mse = evaluate_model(model, test_data, Y_test, hparams["training_type"])

    if bool(hparams['save_model']):
        if hparams['training_type']=='DL':
            model.save(hparams['model_path']+'_{}.h5'.format(mse))
        else:
            with open(hparams['model_path']+'.pkl', 'wb') as file:
                pickle.dump(model, file)
