% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 1. Generate the graph

N = 6;

W = [ 0, 1, 0, 0, 0, 0; ...
      1, 0, 0, 1, 0, 0; ...
      0, 0, 0, 0, 1, 0; ...
      0, 1, 0, 0, 1, 1; ...
      0, 0, 1, 1, 0, 1; ...
      0, 0, 0, 1, 1, 0];

xy = [ 0, 0; ...
       2, 0; ...
       2, 2; ...
       4, 0; ...
       4, 2; ...
       5, 1 ];

% specify a signal
f = (1:N)';

G = GraphSig(W,xy,f,[],'size40 CLim[-1, 1] linewidth2 verbatim{{axis off;}}');


%% 2. Partition the graph

GP = PartitionTreeFiedler(G);


%% 3. Analyze the graph

dmatrix = HGLET_Analysis(G,GP);


%% 4. Specify a basis and visualize it

BS = BasisSpec(1+[1;3;3;2]);
dvec = dmatrix2dvec(dmatrix,GP,BS);
BasisVisual(G,GP,BS,dvec);