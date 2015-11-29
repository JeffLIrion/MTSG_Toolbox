function pm = PartitionPathRatioCut(W)
% Given a weight matrix 'W' for a path graph, partition it using RatioCut
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


%% RatioCut
for row = 1:N-1
    RatioCut = W(row,row+1)*(1/row + 1/(N-row));
    if row == 1
        minRatioCut = RatioCut;
        cutme = 1;
    elseif RatioCut < minRatioCut
        minRatioCut = RatioCut;
        cutme = row;
    end
end

pm = ones(N,1);
pm(cutme+1:N) = -1;


end