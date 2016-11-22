%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         % 
%   Code accompanying                                                     %
%   "Prospective identification of hematopoietic lineage choice by deep   %
%   learning", Nature methods, in Revision                                %
%                                                                         %
%   Authors:                                                              %
%   Felix Buggenthin*, Florian Buettner*, Philipp S Hoppe, Max Endele,    %
%   Manuel Kroiss, Michael Strasser, Michael Schwarzfischer,              % 
%   Dirk Loeffler, Konstantinos D Kokkaliaris, Oliver Hilsenbeck,         %
%   Timm Schroeder†, Fabian J Theis1†, Carsten Marr†                      %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% paths
addpath(genpath('./celldetection'));
%change to the location of the dataset
params.datapath = 'F:\Data\Philipp\FatePrediction_2013\NatMethDemo\';

%% parameters
params.experimentname = 'experiment3';
params.outpath = [params.datapath 'quantification/'];
mkdir(params.outpath)

%% microns per pixel (depends on microscope setup)
params.mperp = 1.0238;


params.doplot = 0;


%segmentation
params.delta = 40;
params.minsize = 20;
params.maxsize = 4000;

params.eccfilter = 0.99;
params.ecccombfilter = 0.7;
params.maxsizecombfilter = 700;
params.minsize = 40;
params.maxsize = 1000;

%cell identification
params.constwindowsize = 30;
params.windowsize = params.constwindowsize;
params.maxwindowsize = 4*params.constwindowsize;
params.dowatershed = 0;
params.growwindow = 10;

params.normalize = 0;

params.patchsize = 41;

params.imagesize = [1388,1040];

params.bfwl = 0;
params.segmentTotal = 1;
%% cell identification
%load tracks
load([params.datapath 'demoTracks.mat'])
unicells = unique(track.cellNr);
%%
for c = unicells
    %compute the cell or load it if it was already computed
    if ~exist([params.outpath '/cell_' num2str(c) '.mat'],'file')
        [Is_c_centered,Iorgs_c_centered,bws_c_centered,cellsizes, cellspeeds, type, label] = identifyCellsFromTracks(track,c,params);
        s = warning('error', 'MATLAB:save:sizeTooBigForMATFile');
        save([params.outpath '/cell_' num2str(c) '.mat'],'Is_c_centered','Iorgs_c_centered','bws_c_centered','cellsizes','cellspeeds', 'label','params','-v7.3');
        disp('Saved.')
    else
        fprintf('Found Cell %i\n',c)
        load([params.outpath '/cell_' num2str(c) '.mat'])
    end
    
    
    
end

