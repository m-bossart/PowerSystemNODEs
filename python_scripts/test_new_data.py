import pandas as pd
import numpy as np
from sklearn.preprocessing import PolynomialFeatures, StandardScaler, MinMaxScaler
from sklearn import linear_model
import matplotlib.pyplot as plt
import tensorflow as tf

from tensorflow import keras
from tensorflow.keras import layers

data = pd.read_excel("IVupdated.xlsx")

def prepare_data(data):
    voltage = list(data["Vm"])
    theta = list(data["Vtheta"])
    time_stamp = list(data["t"].dropna())
    x_train = [np.array([i,j, k]) for i, j, k in zip(voltage, theta, time_stamp)]
    y_train = list(data["Ir"].dropna())
    #count = 1
    #temp_volt = []
    #temp_theta = []
    #for volt, the in zip(voltage, theta):
    #    temp_volt.append(float(volt))
    #    temp_theta.append(float(the))
    #    if count % 100 == 0:
    #        x_train.append([float(sum(temp_volt)/len(temp_volt)), float(sum(temp_theta)/len(temp_theta))])
    #        temp_volt = []
    #        temp_theta = []
    #    count+=1

    #import pdb;pdb.set_trace()
    return np.array(x_train), np.array(y_train), time_stamp

def plot_graphs(ground_truth, predictions, time_stamp):
    fig,ax = plt.subplots()
    ax.plot(time_stamp, ground_truth, color="red", marker="o")
    ax.set_xlabel("Time",fontsize=14)
    ax.set_ylabel("ground_truth",color="red",fontsize=14)

    #ax2=ax.twinx()
# make a plot with different y-axis using second axis object
    ax.plot(time_stamp, predictions ,color="blue",marker="o")
    #ax2.set_ylabel("predictions",color="blue",fontsize=14)
    plt.show()
# save the plot as a file
    fig.savefig('ground truth vs predicted.jpg',
            format='jpeg',
            dpi=100,
            bbox_inches='tight')
    

def lr(x_train, y_train, degree):
    poly = PolynomialFeatures(degree=10)
    x_train = poly.fit_transform(x_train)

    clf = linear_model.LinearRegression()
    clf.fit(x_train, y_train)

    predictions = clf.predict(x_train)
    return predictions

def neural_network(x_train, y_train, train, labels):
    print(x_train.shape)
    labels = labels.reshape(-1,1)
    y_train = y_train.reshape(-1,1)
    scaler_1 = StandardScaler().fit(train)
    #scaler = MinMaxScaler()
    scaler_2 = StandardScaler().fit(labels)
    x_train = scaler_1.transform(x_train)
    y_train = scaler_2.transform(y_train)
    train = scaler_1.transform(train)
    model = tf.keras.Sequential()

    model.add(layers.Dense(128, input_dim=3))
    model.add(layers.Dense(128, activation='tanh'))
    model.add(layers.Dense(128, activation='tanh'))
    model.add(layers.Dense(128, activation='tanh'))
    #model.add(layers.Dense(128, activation='tanh'))
    #model.add(layers.Dense(128, activation='tanh'))
    #model.add(layers.Dense(128, activation='tanh'))
    model.add(layers.Dense(1, activation="linear"))

    model.compile(loss='mean_squared_error',
                optimizer=tf.keras.optimizers.RMSprop(0.001))
    model.fit(x_train, y_train, validation_split=0.0,verbose=0, epochs=100, batch_size=1)
    predictions = model.predict(train)
    #import pdb;pdb.set_trace()
    asd = scaler_2.inverse_transform(predictions, copy=None)
    #print(predictions,asd)
    #return predictions
    return asd

def build_rnn(x_train, y_train, train, labels):
    labels = labels.reshape(-1,1)
    y_train = y_train.reshape(-1,1)
    scaler_1 = StandardScaler().fit(train)
    #scaler = MinMaxScaler()
    scaler_2 = StandardScaler().fit(labels)
    x_train = scaler_1.transform(x_train)
    y_train = scaler_2.transform(y_train)
    train = scaler_1.transform(train)

    x_train = np.expand_dims(x_train, axis=2)
    model = tf.keras.Sequential()
    model.add(layers.Input(shape=(3,1)))
    model.add(layers.SimpleRNN(64, return_sequences=True))
    model.add(layers.SimpleRNN(64, return_sequences=True))
    model.add(layers.SimpleRNN(64, return_sequences=True))
    model.add(layers.SimpleRNN(64, return_sequences=True))
    model.add(layers.SimpleRNN(64, return_sequences=True))
    model.add(layers.SimpleRNN(64, return_sequences=True))
    model.add(layers.SimpleRNN(64))
    model.add(layers.Dense(1, activation="linear"))
    model.compile(loss='mean_absolute_error',
                optimizer=tf.keras.optimizers.Adam(0.0001))
    model.fit(x_train, y_train, validation_split=0.0,verbose=0, epochs=100, batch_size=1, shuffle=False)
    dot_img_file = 'model_1.png'
    tf.keras.utils.plot_model(model, to_file=dot_img_file, show_shapes=True)
    test = np.expand_dims(train, axis=2)
    predictions = model.predict(test)
    asd = scaler_2.inverse_transform(predictions, copy=None)
    #print(predictions, asd)
    #return predictions
    return asd


x_train, y_train, time_stamp = prepare_data(data=data)

y_train = np.array([np.array(i) for i in y_train])
#predictions = neural_network(x_train, y_train)
#import pdb;pdb.set_trace()
counter = 25
number = 20

x_train_trans = x_train[:number]
y_train_trans = y_train[:number]

for _ in range(counter):
    x_train_trans = np.concatenate((x_train_trans, x_train[:number]))
    y_train_trans = np.concatenate((y_train_trans, y_train[:number]))

x_train_trans = np.concatenate((x_train_trans, x_train[number:]))
y_train_trans = np.concatenate((y_train_trans, y_train[number:]))
#import pdb;pdb.set_trace()
#time_stamp = 
predictions = build_rnn(x_train_trans, y_train_trans, x_train, y_train)
#predictions = neural_network(x_train_trans, y_train_trans, x_train, y_train)

plot_graphs(y_train, predictions, time_stamp)
#import pdb;pdb.set_trace()