# HSPC Fate Prediction

Code and data accompanying the publication

**Prospective identification of hematopoietic lineage choice by deep learning**\n
by Felix Buggenthin\*, Florian Buettner\*, Philipp S Hoppe, Max Endele, Manuel Kroiss, Michael Strasser, Michael Schwarzfischer, Dirk Loeffler, Konstantinos D Kokkaliaris, Oliver Hilsenbeck, Timm Schroeder†, Fabian J Theis†, Carsten Marr†*, in Revision.
 
 
## Cell detection (MATLAB)
 Required toolboxes:
 - Image processing toolbox
 - Statistics toolbox
 
 1. Download the example dataset (two positions of experiment 3, full duration, ~10 GB):
 <LINK MISSING>
 2. adjust the path to the dataset in celldetection_metascript.m
 3. Execute celldetection metascript
 
## Predicting cell-specific lineage scores
Required software:
* caffe
* python 2.7
* theano, h5py

Based on the image patches generated using the celldetection script along with the displacemnt feature, our models can be applied to obtain cell-specific predictions of lineage choice. We illustrate the workflow in an ipython notebook that can be viewed [interactively]().  This workflow includes processing of image patches, the extraction of convoluational neural network (CNN)-based patch-specific features as well as the final prediction of lineage choice using a recurent neural network (RNN).
 
## Training
Required software:
* caffe
* python 2.7
* theano, scikit-learn, h5py

In addition to Model training is perfromed in two steps. First, a CNN is trained based on the image patched generated using the celldetection script along with the displacemnt feature.
We provide the caffe model specification for training the model in `CNN_train_test.prototxt` which, along with the solver specifications detaied in `CNN_solver.prototxt`can be used to train the CNN. We provide fully trained models (for all rounds (?)) and solverstates, allowing users to fine-tune models for specific applications. 
After training, the CNNs are used to derive patch-specific features. We provide the extracted features for all experiments as a resource that can be downloaded [here]().

 Next, these CNN-based features are used as input for training an RNN in order to obtain cell-specific lineage scores. 
 RNN training is illustrated in the python script `train_conv.py`. 




