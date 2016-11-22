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
def load_pickle(pickle_file):

	fpickle = open(pickle_file, 'r')
	input_image = pickle.load(fpickle)
	lab = pickle.load(fpickle)
	mov = pickle.load(fpickle)
	cellID = pickle.load(fpickle)
	fpickle.close()

	res = {}
	res['mov'] = mov
	res['im'] = input_image
	res['label'] = lab
	res['cellIDs'] = 'Cell'+cellID[0][0].split('_')[5][4:]
	return res



def load_mat(mat_files):
	res = {}
	res['mov'] = []
	res['im'] = []
	res['label'] = []
	res['cellIDs'] = []


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

		res['mov'].append(mov_test)
		res['im'].append(X)
		res['label'].append(lab[0])
		res['cellIDs'].append(mat_file.split('/')[2].split('.')[0])


	return res


def plotPerfromance(predictions, labels, curve_type='ROC'):
	assert(curve_type in ['ROC', 'PR'])
	
	if curve_type=='ROC':
		auc = metrics.roc_auc_score(labels, predictions) 

		fpr, tpr, _ = metrics.roc_curve(labels, predictions)
		roc_auc = metrics.auc(fpr, tpr)
		plt.plot(fpr, tpr, lw=1, color='navy',
		         label='ROC curve')			
		plt.xlabel('FPR')
		plt.ylabel('TPR')
		plt.title('AUC = {0:0.2f}'.format(auc))

	else:
		precision, recall, _ = metrics.precision_recall_curve(labels,predictions)
		plt.plot(recall, precision, lw=1, color='navy',
		         label='Precision-Recall curve')
		plt.xlabel('Recall')
		plt.ylabel('Precision')


	plt.ylim([0.0, 1.05])
	plt.xlim([0.0, 1.0])
	plt.legend(loc="lower left")
	plt.show()



if __name__ =="__main__":
	film_tr = sys.argv[1]
	num_iter = sys.argv[2]	
	mat_file = '../data/single_cells/cell_1.mat'
	PRETRAINED = '../pretrained/'+film_tr+'/FBalex_BNS_iter_'+str(num_iter)+'.caffemodel'
	MODEL_FILE = 'single_cells/FBalex_deploy.prototxt'
	ims = mat2dict(mat_file)
	res = extract(pretrained=PRETRAINED, model_file=MODEL_FILE, input_image_test=ims['im'], displacement_test=ims['mov'])


