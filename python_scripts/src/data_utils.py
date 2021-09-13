import h5py
import time
import numpy as np

def get_data(filename):
    with h5py.File(filename, 'r') as hdf:
            print("loading Development set (train)...")
            W_dev = np.array(hdf.get('W_dev'))             # W
            X_s_dev = np.array(hdf.get('X_s_dev'))         # X_s
            X_v_dev = np.array(hdf.get('X_v_dev'))         # X_v
            T_dev = np.array(hdf.get('T_dev'))             # T
            Y_dev = np.array(hdf.get('Y_dev'))             # RUL
            A_dev = np.array(hdf.get('A_dev'))             # Auxiliary

            # Test set
            print("loading Test set...")
            W_test = np.array(hdf.get('W_test'))           # W
            X_s_test = np.array(hdf.get('X_s_test'))       # X_s
            X_v_test = np.array(hdf.get('X_v_test'))       # X_v
            T_test = np.array(hdf.get('T_test'))           # T
            Y_test = np.array(hdf.get('Y_test'))           # RUL
            A_test = np.array(hdf.get('A_test'))           # Auxiliary

            # Varnams
            W_var = np.array(hdf.get('W_var'))
            X_s_var = np.array(hdf.get('X_s_var'))
            X_v_var = np.array(hdf.get('X_v_var'))
            T_var = np.array(hdf.get('T_var'))
            A_var = np.array(hdf.get('A_var'))

            # from np.array to list dtype U4/U5
            W_var = list(np.array(W_var, dtype='U20'))
            X_s_var = list(np.array(X_s_var, dtype='U20'))
            X_v_var = list(np.array(X_v_var, dtype='U20'))
            T_var = list(np.array(T_var, dtype='U20'))
            A_var = list(np.array(A_var, dtype='U20'))

    W = np.concatenate((W_dev, W_test), axis=0)
    X_s = np.concatenate((X_s_dev, X_s_test), axis=0)
    X_v = np.concatenate((X_v_dev, X_v_test), axis=0)
    T = np.concatenate((T_dev, T_test), axis=0)
    Y = np.concatenate((Y_dev, Y_test), axis=0)
    A = np.concatenate((A_dev, A_test), axis=0)

    train_data = np.concatenate((W_dev, X_s_dev, X_v_dev, T_dev, A_dev), axis=1)
    test_data = np.concatenate((W_test, X_s_test, X_v_test, T_test, A_test), axis=1)
    all_data = np.concatenate((W, X_s, X_v, T, A), axis=1)
    features_str = sum([i for i in (W_var, X_s_var, X_v_var, T_var, A_var)], [])

    return (train_data, Y_dev), (test_data, Y_test), (all_data, Y), features_str
