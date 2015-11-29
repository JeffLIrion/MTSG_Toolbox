% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~exist('PieceRegular.mat','file')
    fprintf('\n\n\nIn order to generate these results, install WaveLab (<a href="http://statweb.stanford.edu/~wavelab/">http://statweb.stanford.edu/~wavelab/</a>),\nrun WavePath.m, then run this script again.\n\n\n');
    return
end

load('PieceRegular.mat');
G0 = G;
load('PieceRegular_SNR20.mat');

% after the first segment & denoise iteration
GP = PartitionTreePathNCut(G);
[dH,dHrw,dHsym] = HGLET_Path_Analysis_All(G,GP);
[dvec,BS,trans] = HGLET_GHWT_Path_BestBasis_MDL(dH,dHrw,dHsym,0,GP);
[~,G1] = HGLET_GHWT_Path_Synthesis(dvec,GP,BS,trans,G);
Path_BasisVisual(G1,GP,BS,trans);

% segment & denoise and show the result
[G,GP,BS,trans,q,T] = Simultaneous_Segmentation(G);

% info about the run
fprintf('\n\nPiece-Regular:\n');
fprintf('N = %d\n',length(G));
fprintf('Original SNR = 20.00 dB\n');
fprintf('Final SNR = %.2f dB\n',SNR(G0,G));
fprintf('Precision delta = 2^-%d\n',q);
fprintf('Threshold T = %d\n\n',T);