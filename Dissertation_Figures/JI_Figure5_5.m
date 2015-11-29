% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 1. Generate the graph

N = 6;
G = Gpath(N,(1:N)');
[W,xy,f] = ExtractData(G);
W(2,3) = 0.1; W(3,2) = 0.1;
G = GraphSig(W,xy,f);
G = EditPlotSpecs(G,'size250 verbatim{{axis off;}}');


%% 2. Partition the graph

GP = PartitionTreeFiedler(G);


%% 3. Analyze the graph

dmatrix = GHWT_Analysis(G,GP);


%% 4. Specify the yellow basis (fine-to-coarse) and visualize it

BS_yellow = BasisSpec(4-[1;0;0;1;1],[],false);
dvec = dmatrix2dvec(dmatrix,GP,BS_yellow);
BasisVisual(G,GP,BS_yellow,dvec);


%% 5. Specify the green (Haar) basis (fine-to-coarse) and visualize it

BS_green = HaarBasisSpec(GP);
dvec = dmatrix2dvec(dmatrix,GP,BS_green);
BasisVisual(G,GP,BS_green,dvec);