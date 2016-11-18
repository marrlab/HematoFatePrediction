"""
" Multilayer Bidirectional RNN with FastDropout or LSTM units
" Manuel Kroiss <kroissm@in.tum.de>
"""
import os
os.environ["THEANO_FLAGS"] = "mode=FAST_RUN,device=cpu,floatX=float32"

import numpy as np
import theano
import theano.tensor as T
import pickle
import scipy
import time
import matplotlib as mpl
mpl.use('Agg')
from pylab import *
from collections import OrderedDict
import scipy.optimize
from sklearn.metrics import roc_curve, auc, f1_score
import pdb

mode = theano.Mode(linker='cvm', optimizer='fast_compile')
#mode = theano.Mode(linker='py', optimizer='fast_compile')
#mode = 'DEBUG_MODE'

class RNN(object):
		
	def __init__(self, n_input, n_hidden, n_output, architecture):
		self.x = T.tensor3(name='x')
		self.l = T.vector(name='l', dtype='int32')

		self.structure = [n_input, n_hidden, n_output]
		self.params = []
		
		self.rng = np.random.RandomState(1234)
		
		# ff layer
		if architecture == 'ff':
			ha_forw = self.ff_layer(self.x, n_input, n_hidden)
			hb_forw = self.ff_layer(ha_forw, n_hidden, n_hidden)
			h1_forw = self.ff_layer(hb_forw, n_hidden, n_hidden)
		
		# fastdrop layer
		if architecture == 'fastdrop':
			h1_forw = self.fastdrop_recurrent_layer(self.x, n_input, n_hidden)
		
		# lstm layers
		if architecture == 'dblstm':
			h0_forw = self.lstm_layer(self.x, n_input, n_hidden)
			h1_forw = self.lstm_layer(h0_forw, n_hidden, n_hidden)
			h0_back = self.lstm_layer(self.x, n_input, n_hidden, True)
			h1_back = self.lstm_layer(h0_back, n_hidden, n_hidden)

		# last merging layer
		if 'h1_back' in locals():
			self.W_hy = self.make_param(np.random.normal(0, 1, (2 * n_hidden, n_output)) * 0.1)
			self.b_y = self.make_param(np.zeros((n_output,)))
			
			def merge_step(f, b, l):
				h = T.concatenate([f[T.arange(l)], b[T.arange(l-1, -1, -1)]], 1)
				y = T.nnet.sigmoid(T.dot(h, self.W_hy) + self.b_y)
				return T.mean(y, 0)
				
			h_forw_swap = h1_forw.dimshuffle(1, 0, 2)
			h_back_swap = h1_back.dimshuffle(1, 0, 2)
			self.y_mean, _ = theano.map(merge_step, [h_forw_swap, h_back_swap, self.l])
		
		else:
			self.W_hy = self.make_param(np.random.normal(0, 1, (n_hidden, n_output)) * 0.1)
			self.b_y = self.make_param(np.zeros((n_output,)))
			
			def merge_step(f, l):
				h = f[T.arange(l)]
				y = T.nnet.sigmoid(T.dot(h, self.W_hy) + self.b_y)
				#y = fast_dropout_sigmoid_layer(h, self.W_hy, self.b_y, 0.5)
				return T.mean(y, 0)
				
			h_forw_swap = h1_forw.dimshuffle(1, 0, 2)
			self.y_mean, _ = theano.map(merge_step, [h_forw_swap, self.l])
		
		# loss
		self.z = T.matrix(name='z', dtype='int32')
		self.loss = T.mean(T.nnet.binary_crossentropy(self.y_mean * 0.9999 + 0.00005, self.z))
		self.error = T.mean(T.neq(T.round(self.y_mean), self.z))
		self.cost = self.loss

	def prep_dataset(self, data_x, data_y):
		set_x = []
		set_y = []
		
		for i in range(len(data_x)):
			x = data_x[i]
				
			set_x.append(x)
			if len(data_y) > 0:
				y = data_y[i]
				set_y.append(y)

		# length normalize sequences
		lengths = np.array([X.shape[0] for X in set_x], dtype='int32')
		maxlen = np.max(lengths)
		for i in range(len(set_x)):
			set_x[i] = np.resize(set_x[i], (maxlen, set_x[i].shape[1]))

		mat_y = np.array(set_y, dtype=self.z.dtype)
		mat_y = np.reshape(mat_y, (mat_y.shape[0], 1))

		print('OK')
		
		shared_x = theano.shared(np.swapaxes(np.array(set_x, dtype=theano.config.floatX), 0, 1))
		shared_y = theano.shared(mat_y)
		shared_l = theano.shared(lengths)

		return shared_x, shared_l, shared_y

		
	def save(self, filename):
		weights = self.get_theta()
		
		f = open(filename, 'wb')
		pickle.dump([self.structure, weights], f)
		f.close()
		
	@staticmethod
	def load(filename, architecture = 'dblstm'):
		f = open(filename, 'rb')
		[structure, weights] = pickle.load(f)
		f.close()
		
		net = RNN(structure[0], structure[1], structure[2], architecture)		
		net.set_theta(weights)
		
		return net
	
	
	def predict(self, data_x):
		[data_x_shared, data_l_shared, _] = self.prep_dataset(data_x, [])
		predict_model = theano.function([], outputs=self.y_mean, givens={
			self.x: data_x_shared,
			self.l: data_l_shared
		}, mode=mode)
		pp = predict_model()
		return pp
	
	def train_rprop(self, train_x, train_y, test_x, test_y, n_epochs=150, pos_step = 1.2, max_step = 50, neg_step = 0.5, min_step = 10e-6, init_step = 0.0125):
		print('preparing datasets..')
		[train_x_shared, train_l_shared, train_y_shared] = self.prep_dataset(train_x, train_y)
		[test_x_shared, test_l_shared, test_y_shared] = self.prep_dataset(test_x, test_y)
		
		print('compiling..')		
		grads = T.grad(self.cost, self.params)
		updates = OrderedDict()
		
		for param, grad in zip(self.params, grads):
			last_grad = theano.shared(param.get_value()*0.)
			last_delta = theano.shared(param.get_value()*0.+init_step)

			change = T.sgn(grad * last_grad)
			onward = T.gt(change, 0)
			turned = T.lt(change, 0)
			delta = T.switch(onward,
				T.minimum(last_delta * pos_step, max_step),
				T.switch(turned,
					T.maximum(last_delta * neg_step, min_step),
					last_delta
			))

			updates[param] = param - T.sgn(grad) * delta
			updates[last_grad] = T.switch(turned, 0, grad)
			updates[last_delta] = delta

		train_model = theano.function([], [self.cost, self.error], givens={
			self.x: train_x_shared,
			self.l: train_l_shared,
			self.z: train_y_shared
		}, updates=updates, mode=mode)
		
		test_model = theano.function([], [self.loss, self.error, self.y_mean], givens={
			self.x: test_x_shared,
			self.l: test_l_shared,
			self.z: test_y_shared
		}, mode=mode)

		best_cost = float('Inf')
		best_params = None
		monitor = []
		monitor_param = []

		print('training..')
		for epoch in range(n_epochs):
			start = time.time()
			
			[cst, err] = train_model()
			print('train	epoch={:03d}   loss={:.5f}   error={:.5f}'.format(epoch+1, float(cst), float(err)))
			
			[cst2, err2, ym2] = test_model()
			p = ym2.flatten()
			s = np.array(test_y).flatten()
			f1 = f1_score(s, (np.array(p)>.5)*1,pos_label = None, average='macro')
			fpr, tpr, thresholds = roc_curve(s, p)
			roc_auc = auc(fpr, tpr)
			print('	 valid	epoch={:03d}   loss={:.5f}   error={:.5f}	AUC={:.5f} F1={:.5f}'.format(epoch+1, float(cst2), float(err2), float(roc_auc), float(f1)))
			
			print('	 elapsed={:.5f}'.format(time.time()-start))
			monitor.append([cst2, err, err2, roc_auc, f1])
			monitor_param.append(self.get_theta())

			if cst2 < best_cost:
				best_cost = cst2
				best_params = self.get_theta()
				weights = self.get_theta()
			
			#fn = 'model1214Conv_'+str(epoch)+'.plk'
			#f = open(fn, 'wb')
			#pickle.dump([self.structure, weights], f)
			#f.close()	
			self.set_theta(best_params)
		
		plot(range(len(monitor)), np.array(monitor))
		show()
		res = {}
		res['perf'] = monitor
		res['param'] = monitor_param
		return res

	def get_theta(self):
		return np.concatenate([ p.get_value().ravel() for p in self.params ])
	
	def set_theta(self, theta):
		offset = 0
		for param in self.params:
			shape = param.shape.eval()
			len = np.prod(shape)
			param.set_value(theta[offset:offset+len].reshape(shape).astype('float32'))
			offset += len
			
	  
	def make_param(self, init):
		param = theano.shared(np.array(init, dtype=theano.config.floatX))
		self.params.append(param)
		return param
		
	def ff_layer(self, x, n_in, n_out, irange=0.05):
		W_xh = self.make_param(np.random.uniform(-irange, irange, (n_in, n_out)))
		b_h = self.make_param(np.zeros((n_out,)))
		
		output = T.tanh(T.dot(x, W_xh) + b_h)
		return output

	def fastdrop_recurrent_layer(self, x, n_in, n_out, go_backwards=False):
		W_hh_init = np.random.normal(0, 1, (n_out, n_out)) * 0.1
		for i in range(n_out):
			W_hh_init[i][np.random.choice(n_out, n_out-16)] = 0
		pWhh = np.max(abs(scipy.linalg.eigvals(W_hh_init)))
		W_hh_init = W_hh_init / pWhh * 1.1
		
		W_xh = self.make_param(np.random.normal(0, 1, (n_in, n_out)) * 0.001)
		W_hh = self.make_param(W_hh_init)
		h_0 = self.make_param(np.zeros((n_out,)))
		b_h = self.make_param(np.zeros((n_out,)))
			
		def step(x_t, h_tm1):
			c_t = T.concatenate([x_t, h_tm1], 1)
			W_ch = T.concatenate([W_xh, W_hh])
			h_t = fast_dropout_tanh_layer(c_t, W_ch, b_h, 0.3)
			return h_t
			
		output, _ = theano.scan(step, sequences=x, outputs_info=[T.alloc(h_0, x.shape[1], n_out)], go_backwards=go_backwards)
		return output		
		
	def lstm_layer(self, x, n_in, n_out, go_backwards=False):
		W_xh = self.make_param(np.random.uniform(-0.1, 0.1, (n_in, n_out*4)))
		W_hh = self.make_param(np.random.uniform(-0.1, 0.1, (n_out, n_out*4)))
		b_h = self.make_param(np.zeros((n_out*4, )))
		W_ci = self.make_param(np.random.uniform(-0.1, 0.1, (n_out, )))
		W_cf = self.make_param(np.random.uniform(-0.1, 0.1, (n_out, )))
		W_co = self.make_param(np.random.uniform(-0.1, 0.1, (n_out, )))
		c_0 = self.make_param(np.zeros((n_out, )))
		h_0 = self.make_param(np.zeros((n_out, )))

		def step(x_t, c_tm1, h_tm1):
			pre = T.dot(x_t, W_xh) + T.dot(h_tm1, W_hh) + b_h
			i_t = T.nnet.sigmoid(pre[:, :n_out] + c_tm1 * W_ci)
			f_t = T.nnet.sigmoid(pre[:, n_out:2*n_out] + c_tm1 * W_cf)
			c_t = f_t * c_tm1 + i_t * T.tanh(pre[:, 2*n_out:3*n_out])
			o_t = T.nnet.sigmoid(pre[:, 3*n_out:4*n_out] + c_t * W_co)
			h_t = o_t * T.tanh(c_t)
			return c_t, h_t
			
		[_, output], _ = theano.scan(step, sequences=x, outputs_info=[T.alloc(c_0, x.shape[1], n_out), T.alloc(h_0, x.shape[1], n_out)], go_backwards=go_backwards)
		return output
	
	
def fast_dropout_tanh_layer(x, W, b, p, alpha=1):
	""" p is the probablity of keeping a unit
	"""
	mu2_x = T.square(T.mean(x, 0))
	si2_x = T.mean(T.square(x), 0) - mu2_x
			
	s2 = T.dot(alpha * p * (1-p) * mu2_x + p * si2_x, T.square(W))
	mu = p * T.dot(x, W) + b
	pre = (2 * mu) / T.sqrt(1 + 0.125 * np.pi * (4 * s2))
	h = 2 * T.nnet.sigmoid(pre) - 1
		
	return h
