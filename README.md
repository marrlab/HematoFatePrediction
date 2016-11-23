# Hemato Fate Prediction

Code and data accompanying 

**Prospective identification of hematopoietic lineage choice by deep learning**

by Felix Buggenthin\*, Florian Buettner\*, Philipp S Hoppe, Max Endele, Manuel Kroiss, Michael Strasser, Michael Schwarzfischer, Dirk Loeffler, Konstantinos D Kokkaliaris, Oliver Hilsenbeck, Timm Schroeder†, Fabian J Theis†, Carsten Marr† 

Download the data at https://hmgubox.helmholtz-muenchen.de:8001/d/ccbfb5f1ac/
 
 
## Cell detection
Required software:
* MATLAB (R2014a)
* MATLAB Image processing toolbox
* MATLAB Statistics toolbox

Steps:
 1. Download the example dataset (two positions of experiment 3, full duration, ~10 GB):
 <LINK MISSING>
 2. adjust the path to the dataset in celldetection_metascript.m
 3. Execute celldetection metascript
 
## Predicting cell-specific lineage scores
Required software:
* caffe ([this fork](https://github.com/flophys/caffe) allowing for prediction with concatenation layer) 
* python 2.7
* theano>=0.8.2, scikit-learn>=0.18.1, h5py>=2.6.0 

Based on the image patches generated using the celldetection script along with the displacemnt feature, our models can be applied to obtain cell-specific predictions of lineage choice. We illustrate the workflow in an ipython notebook that can be viewed [interactively](http://nbviewer.ipython.org/github/QSCD/HematoFatePrediction/blob/master/cellprediction/Predict_cell_fates.ipynb).  This workflow includes processing of image patches, the extraction of convoluational neural network (CNN)-based patch-specific features as well as the final prediction of lineage choice using a recurent neural network (RNN).
 
## Training
Required software:
* caffe ([this fork](https://github.com/flophys/caffe) for prediction with concatenation layer) 
* python 2.7
* theano>=0.8.2, scikit-learn>=0.18.1, h5py>=2.6.0

To install caffe, please follow these [installation instructions](http://caffe.berkeleyvision.org/installation.html) for your OS. We highly recommend using the [Anaconda framework](https://docs.continuum.io).  


In addition to Model training is perfromed in two steps. First, a CNN is trained based on the image patched generated using the celldetection script along with the displacemnt feature.
We provide the caffe model specification for training the model in `CNN_train_test.prototxt` which, along with the solver specifications detaied in `CNN_solver.prototxt`can be used to train the CNN. We further provide a fully trained model and solverstate, allowing users to fine-tune models for specific applications. After training, the CNN is used to derive patch-specific features. We provide the extracted features for all experiments as a resource that can be downloaded [here]().

 Next, these CNN-based features are used as input for training an RNN in order to obtain cell-specific lineage scores. 
 RNN training is illustrated in the python script [`train_conv.py`](https://github.com/QSCD/HematoFatePrediction/blob/master/cellprediction/py/train_conv.py). 




