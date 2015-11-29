% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% the original "Blocks"
load('Blocks.mat');
[~,xy] = ExtractData(G);
xmin = min(xy)-1;
xmax = max(xy)+1;
G = EditPlotSpecs(G,sprintf('notitle linewidth3 verbatim{{xlim([%f,%f])}}',xmin,xmax));
GraphSig_Plot(G);
G0 = G;

% the noisy "Blocks"
load('Blocks_Noisy.mat');
G = EditPlotSpecs(G,sprintf('notitle linewidth3 verbatim{{xlim([%f,%f])}}',xmin,xmax));
GraphSig_Plot(G);

% info about the signals
fprintf('\n\nBlocks:\n');
fprintf('N = %d\n',length(G));
fprintf('SNR = %.2f dB\n\n',SNR(G0,G));