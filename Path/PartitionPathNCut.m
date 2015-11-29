function pm = PartitionPathNCut(W)
% Given a weight matrix 'W' for a path graph, partition it using Ncut
%
% Inputs
%   W           the edge weight matrix
%
% Outputs
%   pm          a vector of 1's and -1's
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% Preliminaries

N = length(W);


%% handle the case when there are 2 nodes
if N == 2
    pm = [1;-1];
    return
end


    
%% Ncut
D = sum(W);
for row = 1:N-1
    Ncut = W(row,row+1)*(1/sum(D(1:row)) + 1/sum(D(row+1:N)));
    if row == 1
        minNcut = Ncut;
        cutme = 1;
    elseif Ncut < minNcut
        minNcut = Ncut;
        cutme = row;
    end
end

pm = ones(N,1);
pm(cutme+1:N) = -1;


end