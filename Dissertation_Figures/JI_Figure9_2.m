% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~exist('Msignal.mat','file')
    fprintf('\n\n\nIn order to generate these results, install WaveLab (<a href="http://statweb.stanford.edu/~wavelab/">http://statweb.stanford.edu/~wavelab/</a>),\nrun WavePath.m, then run this script again.\n\n\n');
    return
end

load('Msignal.mat');

% plot the original signal
[~,xy] = ExtractData(G);
xmin = min(xy)-1;
xmax = max(xy)+1;
G = EditPlotSpecs(G,sprintf('notitle linewidth3 verbatim{{xlim([%f,%f])}}',xmin,xmax));
GraphSig_Plot(G);

% segment & denoise and show the result
[G,GP,BS,trans,q,T] = Simultaneous_Segmentation(G);

% info about the run
fprintf('\n\nMsignal:\n');
fprintf('N = %d\n',length(G));
fprintf('Precision delta = 2^-%d\n',q);
fprintf('Threshold T = %d\n\n',T);