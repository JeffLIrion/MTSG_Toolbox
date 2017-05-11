% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminaries
load('Toronto.mat','G');
N = length(G);

% partition the graph
GP = PartitionTreeFiedler(G,1);


%% 1. Analyze the signal

dmatrixG = GHWT_Analysis(G,GP);

% 3) GHWT best basis
[dvec3,BS3,~,p3] = HGLET_GHWT_BestBasis_minrelerror(0,0,0,dmatrixG,GP,G);
r3 = orth2relerror(dvec3);


%% Figure 5
G = EditPlotSpecs(G,'sortnodes size5 verbatim{{axis equal;}}');
BasisVisual(G,GP,BS3,dvec3);