% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



load('Blocks.mat');
G0 = G;
load('Blocks_Noisy.mat');

% specify the "true" segmentation
ind = (1:2048)';
rs = zeros(2049,3);
rs(1,:) = 1;
rs(2,1) = 2049;
rs(:,3) = (1:2049)';
rs(1:13,2) = [1; 206; 267; 308; 472; 513; 820; 902; 1332; 1557; 1598; 1660; 2049];
GP = GraphPart(ind,rs);
levlengths = diff(rs(1:13,2));
levlist = 2*ones(length(levlengths),1);
BS = BasisSpec(levlist,diff(rs));

% display the noisy signal and its starting segmentation
trans = false(length(levlist),2);
Path_BasisVisual(G,GP,BS,trans,round(length(G)/50));
ylims = ylim;

% segment & denoise and show the result
[G,GP,BS,trans,q,T] = Simultaneous_Segmentation(G,round(length(G)/50),GP,BS);
ylim(ylims);

% info about the run
fprintf('\n\nBlocks:\n');
fprintf('N = %d\n',length(G));
fprintf('Original SNR = 11.95 dB\n');
fprintf('Final SNR = %.2f dB\n',SNR(G0,G));
fprintf('Precision delta = 2^-%d\n',q);
fprintf('Threshold T = %d\n\n',T);