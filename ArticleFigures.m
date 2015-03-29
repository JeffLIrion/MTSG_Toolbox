% ArticleFigures
% Generate selected figures from "Hierarchical Graph Laplacian Eigen
% Transforms" and "The Generalized Haar-Walsh Transform"
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% Figures from "Hierarchical Graph Laplacian Eigen Transforms"

%%% FIGURE 1

load('MN_MutGauss.mat');

% convert to inverse Euclidean weights
G = Adj2InvEuc(G);

% partition the graph
GP = PartitionTreeFiedler(G,1);

% extract the regionstarts data
[~,rs] = ExtractData(GP);

% Fig 1a
H1a = HGLET_BasisVector(G,GP,rs(2,2)+4,2,1);

% Fig 1b
H1b = HGLET_BasisVector(G,GP,rs(3,3)+1,3,1);

% Fig 1c
H1c = HGLET_BasisVector(G,GP,rs(6,4)+2,4,1);

% Fig 1d
H1d = HGLET_BasisVector(G,GP,rs(6,4)+9,4,1);

% NOTE: Due to updates in the code since that paper was published, the
% signs of the figures generated herein for Figs. 1c and 1d are opposite of
% those in the article.  



%% Figures from "The Generalized Haar-Walsh Transforms"

%%% FIGURE 1

% Fig 1a
G1a = GHWT_BasisVector(G,GP,2,1);

% Fig 1b
G1b = GHWT_BasisVector(G,GP,8,1);

% Fig 1c
G1c = GHWT_BasisVector(G,GP,3,3);

% Fig 1d
G1d = GHWT_BasisVector(G,GP,rs(8,4)+3,4);


%%% FIGURE 4

% add noise to the mutilated Gaussian
GN = AddNoise(G,5,'gaussian');

% perform the GHWT transform on the noisy signal
dmatrix = GHWT_Analysis(GN,GP);

% find the GHWT best basis
[dvec,BS] = GHWT_BestBasis(dmatrix,GP,0.1);

% threshold the HGLET coefficients
[dvecT,kept] = dvec_Threshold(dvec,'s',0.11,GP,BS);

% reconstruct
[~,G4c] = GHWT_Synthesis(dvecT,GP,BS,GN);

% print the SNR of the reconstructed signal
fprintf('\n\nThe GHWT reconstructed signal has SNR %.2f dB (%d%% coefficients kept)\n\n',SNR(G,G4c),round(100*kept/length(dvecT)));

% take the absolute value of the original minus the denoised signal
G4d = abs(GN-G4c);

% find the min and max values for the colormap
cmin = min([min(G), min(GN), min(G4c), min(G4d)]);
cmax = max([max(G), max(GN), max(G4c), max(G4d)]);

% set the min and max values for the colormap
G4a = EditPlotSpecs(G,sprintf('notitle CLim[%f,%f]',cmin,cmax));
G4b = EditPlotSpecs(GN,sprintf('notitle CLim[%f,%f]',cmin,cmax));
G4c = EditPlotSpecs(G4c,sprintf('notitle CLim[%f,%f]',cmin,cmax));
G4d = EditPlotSpecs(G4d,sprintf('notitle CLim[%f,%f]',cmin,cmax));

% NOTE: These figures and results will differ from those in the article
% because the noise is generated randomly and thus not the same as in the
% experiment that produced the article figures.  

% clear unnecessary variables
clear rs dmatrix dvec dvecT kept cmin cmax