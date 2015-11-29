% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% load the graph
load('Dendrite.mat','G');

% partition the graph (just a formality)
GP = PartitionTreeFiedler(G);

% compute the eigenvector
G1142 = HGLET_BasisVector(G,GP,1143,1);
G1142 = EditPlotSpecs(G1142,'symm notitle size20 sortnodes verbatim{{axis equal; axis off;}}');

% display the figure
GraphSig_Plot(G1142);