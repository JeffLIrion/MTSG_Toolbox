% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



G = Gpath(6);
GP = PartitionTreeFiedler(G);

% coarse-to-fine dictionary
GHWT_Dictionary_c2f(G,GP);

% fine-to-coarse dictionary
GHWT_Dictionary_f2c(G,GP);