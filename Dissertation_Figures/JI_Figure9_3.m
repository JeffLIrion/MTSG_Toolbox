% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~exist('PieceRegular.mat','file')
    fprintf('\n\n\nIn order to generate these results, install WaveLab (<a href="http://statweb.stanford.edu/~wavelab/">http://statweb.stanford.edu/~wavelab/</a>),\nrun WavePath.m, then run this script again.\n\n\n');
    return
end

% the original "Piece-Regular"
load('PieceRegular.mat');
[~,xy] = ExtractData(G);
xmin = min(xy)-1;
xmax = max(xy)+1;
G = EditPlotSpecs(G,sprintf('notitle linewidth3 verbatim{{xlim([%f,%f])}}',xmin,xmax));
GraphSig_Plot(G);
G0 = G;

% the noisy "Piece-Regular"
load('PieceRegular_SNR20.mat');
G = EditPlotSpecs(G,sprintf('notitle linewidth3 verbatim{{xlim([%f,%f])}}',xmin,xmax));
GraphSig_Plot(G);

% info about the signals
fprintf('\n\nPiece-Regular:\n');
fprintf('N = %d\n',length(G));
fprintf('SNR = %.2f dB\n\n',SNR(G0,G));