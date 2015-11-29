function GP = PartitionTreePathRatioCut(G)
% Generate a partition tree for a path graph using RatioCut
%
% Input
%   G       a GraphSig object
%
% Output
%   GP      a GraphPart object
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminary stuff

% constants
W = ExtractData(G);
N = length(W);
jmax = max(3*floor(log2(N)),4);

% TRACKING VARIABLES

% the way in which the nodes are indexed on each level
sslac = ind_class(N);
ind = zeros(N,jmax,sslac);
    
ind(:,1) = (1:N)';

% stands for regionstarts, meaning that the index in 'ind' of the first
% point in region number i is rs(i)
rs = zeros(N+1,jmax,sslac);
rs(1,:) = 1;
rs(2,1) = N+1;



%% 1. Partition the graph to yield rs, ind, and jmax

j = 1;
regioncount = 0;

% use RatioCut to recursively bipartition
while regioncount < N
    % the number of regions on level j
    regioncount = nnz(rs(:,j))-1;

    % add a column for level j+1, if necessary
    if j == jmax
        rs(1,j+1) = 1;
        jmax = jmax+1;
    end

    for r = 1:regioncount
        % regions with 2 or more nodes
        if rs(r,j) ~= rs(r+1,j)-1
            rs1 = rs(r,j);
            rs2 = rs(r+1,j)-1;
            indrs = ind(rs(r,j):rs(r+1,j)-1,j);

            % partition the current region
            pm = PartitionPathRatioCut(W(indrs,indrs));
            r1 = sum(pm > 0);% # points in subregion 1
            r2 = sum(pm < 0);% # points in subregion 2

            % update the indexing
            indr  = zeros(r1+r2,1);
            indr(1:r1)       = indrs(pm > 0);
            indr(r1+1:r1+r2) = indrs(pm < 0);
            ind(rs1:rs2,j+1) = indr;
            clear indr

            % update the region tracking
            rr = nnz(rs(:,j+1));
            rs(rr+1,j+1) = rs1+r1;
            rs(rr+2,j+1) = rs2+1;

        % regions with 1 node
        else
            rr = nnz(rs(:,j+1));
            rs(rr+1,j+1) = rs(r+1,j);
            ind(rs(rr,j+1),j+1) = ind(rs(r,j),j);                    
        end
    end

    % take it to the next level
    j = j+1;
end


% get rid of excess columns in rs
rs(:,j:end) = [];

% get rid of all but the last column of ind
ind = ind(:,j-1);

% create a GraphPart object
GP = GraphPart(ind,rs,[],[],[],[],[],'RatioCut');


end