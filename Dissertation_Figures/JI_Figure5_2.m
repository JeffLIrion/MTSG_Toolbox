% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 1. Generate the graph

N = 6;
G = Gpath(N);
[W,xy,f] = ExtractData(G);
W(2,3) = 0.1; W(3,2) = 0.1;
G = GraphSig(W,xy,f);


%% 2. Partition the graph

GP = PartitionTreeFiedler(G);


%% 3. Display the basis vectors

GHWT_Dictionary_f2c(G,GP);


%% 4. Count the number of choosable bases & display the partitioning tree

[c2f,f2c] = HGLET_GHWT_BasisCount(GP);
fprintf('\n\nThere are:\n\n');
fprintf('* %2d choosable orthonormal bases from the coarse-to-fine dictionary\n',c2f);
fprintf('* %2d choosable orthonormal bases from the fine-to-coarse dictionary\n     (not contained in the coarse-to-fine dictionary)\n',f2c);
fprintf('* %2d total choosable orthonormal bases\n\n',c2f+f2c);

PartitionTreeDisplay(GP,1);
title(sprintf('Fine-To-Coarse Dictionary Structure (%d choosable bases)',f2c));