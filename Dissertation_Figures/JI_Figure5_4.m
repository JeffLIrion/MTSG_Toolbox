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


%% 4. Specify a basis (coarse-to-fine) and visualize it

BS = BasisSpec(1+[2;2;1]);
dvec = dmatrix2dvec(dmatrix,GP,BS);
BasisVisual(G,GP,BS,dvec);