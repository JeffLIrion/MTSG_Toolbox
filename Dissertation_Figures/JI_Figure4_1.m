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

f = zeros(N,1);

G = GraphSig(W,xy,f,[],'size40 CLim[-1, 1] linewidth2');


%% 2. Partition the graph

GP = PartitionTreeFiedler(G);


%% 3. Display the basis vectors

HGLET_Dictionary(G,GP);


%% 4. Count the number of chooseable bases & display the partitioning tree

count = HGLET_GHWT_BasisCount(GP);
fprintf('\n\nThere are:\n\n');
fprintf('* %2d choosable orthonormal bases\n\n',count);

PartitionTreeDisplay(GP);
title(sprintf('Recursive Partitioning Tree <==> Dictionary Structure (%d choosable bases)',count));