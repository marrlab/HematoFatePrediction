# Hemato Fate Prediction

Code and data accompanying 

**Prospective identification of hematopoietic lineage choice by deep learning**

by Felix Buggenthin\*, Florian Buettner\*, Philipp S Hoppe, Max Endele, Manuel Kroiss, Michael Strasser, Michael Schwarzfischer, Dirk Loeffler, Konstantinos D Kokkaliaris, Oliver Hilsenbeck, Timm Schroeder†, Fabian J Theis†, Carsten Marr† 

published in Nature Methods in 2017 (dOi:10.1038/nmeth.4182)

Download the required data from https://drive.google.com/file/d/1j10HeL87CIkdvHzC-98IUt-IUh9XHHTb/view?usp=sharing
 
 
## Cell detection
Required software:
* MATLAB (R2014a)
* MATLAB Image processing toolbox
* MATLAB Statistics toolbox

Steps:
 1. Download the dataset Rawdata_buggenthin_buettner_naturemethods2016 (two exemplary positions of experiment 3, ~10 GB) from the link above
 2. Adjust the path to the dataset in celldetection_metascript.m in our repository
 3. Execute celldetection_metascript.m
 
## Cell prediction
Required software:
* caffe ([this fork](https://github.com/flophys/caffe) allowing for prediction with concatenation layer) 
* python 2.7
* theano>=0.8.2, scikit-learn>=0.18.1, h5py>=2.6.0 

### Predicting lineage scores
Based on the image patches generated using the celldetection_metascript.m along with the displacemnt feature, our models can be applied to obtain cell-specific predictions of lineage choice. We illustrate the workflow in an ipython notebook that can be viewed [interactively](http://nbviewer.ipython.org/github/QSCD/HematoFatePrediction/blob/master/cellprediction/Predict_cell_fates.ipynb).  This workflow includes processing of image patches, the extraction of convoluational neural network (CNN)-based patch-specific features as well as the final prediction of cell-specific lineage scores using a recurent neural network (RNN).
 
### Training the networks
Required software:
* caffe ([this fork](https://github.com/flophys/caffe) for prediction with concatenation layer) 
* python 2.7
* theano>=0.8.2, scikit-learn>=0.18.1, h5py>=2.6.0

To install caffe, please follow these [installation instructions](http://caffe.berkeleyvision.org/installation.html) for your OS. We highly recommend using the [Anaconda framework](https://docs.continuum.io).  


Model training is performed in two steps. First, a CNN is trained based on the image patches generated using the celldetection_metascript.m along with the displacemnt feature.
We provide the caffe model specification for training the model in `CNN_train_test.prototxt` which, along with the solver specifications detaied in `CNN_solver.prototxt` can be used to train the CNN. We further provide a fully trained model and solverstate, allowing users to fine-tune models for specific applications. After training, the CNN is used to derive patch-specific features.

 Next, these CNN-based features are used as input for training an RNN in order to obtain cell-specific lineage scores. 
 RNN training is illustrated in the python script [`train_conv.py`](https://github.com/QSCD/HematoFatePrediction/blob/master/cellprediction/py/train_conv.py). 




