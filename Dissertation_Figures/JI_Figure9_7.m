% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



load('Blocks.mat');
G0 = G;
load('Blocks_Noisy.mat');

% segment & denoise and show the result
[G,GP,BS,trans,q,T] = Simultaneous_Segmentation(G);

% show the result without absorbing small regions
Path_BasisVisual(G,GP,BS,trans,0);

% info about the run
fprintf('\n\nBlocks:\n');
fprintf('N = %d\n',length(G));
fprintf('Original SNR = 11.95 dB\n');
fprintf('Final SNR = %.2f dB\n',SNR(G0,G));
fprintf('Precision delta = 2^-%d\n',q);
fprintf('Threshold T = %d\n\n',T);