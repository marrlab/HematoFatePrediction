from rnn import RNN
import h5py
import numpy as np
import scipy as SP
import scipy.io
import os as os
import cPickle as pickle
def make_sets(X, Y, ratio):
    np.random.seed(1)
    indices = np.random.permutation(len(X))
    cut = len(indices)*ratio
    train, valid = indices[:cut], indices[cut:]
    train_x = [X[i] for i in train]
    train_y = [Y[i] for i in train]
    valid_x = [X[i] for i in valid]
    valid_y = [Y[i] for i in valid]
    return train_x, train_y, valid_x, valid_y
    

film_tr = '120602PH5'
film_tr2 = '130218PH8'


ftr = '../training_data/feat_'+film_tr+'_'+film_tr2+'_TRAIN.pickle'

#training data
file_tr = open(ftr,'r')
res = pickle.load(file_tr)
X = res['feats']
Y = res['lab']
file_tr.close()


## generate validation file
train_x = X 
train_y = Y 
train_x, train_y, valid_x, valid_y = make_sets(train_x, train_y, 0.85)

## train model
print("Train Balance: " + str(len(train_x)) + "*" + str(np.mean(train_y)))
print("Valid Balance: " + str(len(valid_x)) + "*" + str(np.mean(valid_y)))
idim = len(train_x[0][0])
odim = max(train_y)+1

## train model
ndim = len(train_x[0][0])
print("Dimensions: " + str(ndim))

auc_list = list()
f1_list = list()
theta_list = list()


model = RNN(n_input=ndim, n_hidden=20, n_output=1, architecture='dblstm')
res = model.train_rprop(train_x, train_y, valid_x, valid_y, n_epochs=15)
fname = '../models/rnn_models/trained_modelRNN'+film_tr+'test' + '.pkl'
model.save(fname )


