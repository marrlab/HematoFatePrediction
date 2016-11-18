# coding: utf-8
# Classifying with LeNet
import numpy as np
import glob
import pdb
import matplotlib as mpl
#if 1:
import h5py
import scipy as SP
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn import metrics
import random

caffe_root = '../'  # this file is expected to be in {caffe_root}/examples
import sys
sys.path.insert(0, caffe_root + 'python')
import pickle
import caffe


def extract(pretrained, model_file,input_image_test,displacement_test):
#PRETRAINED = '../snapshots/'+film_tr+'A/FBalexNM2_BNS_iter_'+str(num_iter)+'.caffemodel'

	net = caffe.Classifier(model_file, pretrained)

	net.transformer.set_input_scale('data0', 0.00390625)	

	preds_all = list()
	numFeats = 51

	numCells = len(input_image_test)
	featMat = np.zeros((numCells,numFeats))
	pred = np.zeros((numCells,1))
	for num_cell in range(numCells):
		pred[num_cell]= net.predictFB([input_image_test[num_cell]],np.reshape(displacement_test[num_cell],(1,1)), oversample=False)[:,1]#.mean()
	#pred_[num_cell]= net.predict([input_image_test[numim][num_cell]], oversample=False)[:,1]#.mean()
		featMat[num_cell,:50] = net.blobs['fc7_BN'].data.ravel()
		featMat[num_cell,50] = displacement_test[num_cell]#net.blobs['hdfdata'].data.ravel()[0] 

	res = dict()
	res['feats'] = featMat
	res['pred_all'] = pred
	return res

#def extractFeatures(pretrained=PRETRAINED, model_file=MODEL_FILE. im_list = im_list)
#	feat	


def mat2dict(mat_files):
	res_list = []
	for mat_file in mat_files:
		mat = h5py.File(mat_file, 'r')
		matX = mat['Is_c_centered'].value
		matX2 = mat['cellspeeds'][:][0]
		matL = mat['type'][:][()]


		X = []
		X2 = []
		lab = []
		for i in range(len(matX2)):
			if not SP.isnan(matX2[i]):
				oldx = mat[matX[i][0]].value.T
				newx = np.array(oldx, dtype=np.float32)
				newx=newx-newx.mean()
				newx/=newx.std()
				newx = SP.misc.fromimage(SP.misc.toimage(newx)).astype('float') #convert to greyscale image
				newx = SP.array(newx[:,:,SP.newaxis])
				X.append(newx)

				#oldx2 = matX2[i][0,]
				oldx2 =matX2[i]
				newx2 = np.array(oldx2, dtype=np.float32)
				X2.append(newx2) 

				oldlab =matL
				newlab = np.array(oldlab, dtype=np.float32)
				lab.append(newlab) 				
			
		mov_test = np.array(X2)
		mov_test-=mov_test.mean(0)
		mov_test/=mov_test.std(0)
		np.reshape(mov_test, (len(X),1)).shape	

		res = {}
		res['mov'] = mov_test
		res['im'] = X
		res['label'] = lab

		res_list.append(res)
	return res_list

if __name__ =="__main__":
	film_tr = sys.argv[1]
	num_iter = sys.argv[2]	
	mat_file = '../data/single_cells/cell_1.mat'
	PRETRAINED = '../pretrained/'+film_tr+'/FBalex_BNS_iter_'+str(num_iter)+'.caffemodel'
	MODEL_FILE = 'single_cells/FBalex_deploy.prototxt'
	ims = mat2dict(mat_file)
	res = extract(pretrained=PRETRAINED, model_file=MODEL_FILE, input_image_test=ims['im'], displacement_test=ims['mov'])


