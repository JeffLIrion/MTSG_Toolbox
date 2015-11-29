% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



load('Dendrite.mat','G');
G = EditPlotSpecs(G,'notitle sortnodes size10 verbatim{{axis equal; axis off;}}');
GraphSig_Plot(G);