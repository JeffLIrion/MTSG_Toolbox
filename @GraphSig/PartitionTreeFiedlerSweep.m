function GP = PartitionTreeFiedlerSweep(G,~)
% Generate a partition tree for a graph by sweeping through the entries of
% the Fiedler vectors of either L (the unnormalized Laplacian) or L_rw (the
% random-walk normalized Laplacian).  Whereas PartitionTreeFiedler cuts the
% eigenvector entries based on whether they are <0 or >0, this function
% (1) sweeps through the entries of the Fiedler vector and finds the 
% largest difference between successive sorted eigenvector entries, (2) 
% sets "cutoff" to be the midpoint of this interval, and (3) partitions the
% eigenvector entries based on whether they are <cutoff or >cutoff.  
%
% Input
%   G       a GraphSig object
%   ~       if a 2nd input is given, use L_rw for partitioning
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
N = G.length;
jmax = max(3*floor(log2(N)),4);

% TRACKING VARIABLES

% the way in which the nodes are indexed on each level
sslac = ind_class(N);
ind = zeros(N,1,sslac);
ind = (1:N)';

% stands for regionstarts, meaning that the index in 'ind' of the first
% point in region number i is rs(i)
rs = zeros(N+1,jmax,sslac);
rs(1,:) = 1;
rs(2,1) = N+1;



%% 1. Partition the graph to yield rs, ind, and jmax

j = 1;
regioncount = 0;

% L (unnormalized)
if nargin == 1
    while regioncount < N
        % the number of regions on level j
        regioncount = nnz(rs(:,j))-1;
        
        % add a column for level j+1, if necessary
        if j == jmax
            rs(1,j+1) = 1;
            jmax = jmax+1;
        end
        
        % for tracking the child regions
        rr = 1;
        
        % cycle through the parent regions
        for r = 1:regioncount
            % the start of the parent region and 1 node after the end of
            % the parent region
            rs1 = rs(r,j);
            rs2 = rs(r+1,j);
            
            % the number of nodes in the parent region
            n = rs2-rs1;
            
            % regions with 2 or more nodes
            if n > 1
                indrs = ind(rs1:rs2-1);

                % partition the current region
                pm = PartitionFiedlerSweep(G.W(indrs,indrs));
                
                % determine the number of points in child region 1
                n1 = sum(pm > 0);

                % switch regions 1 and 2, if necessary, based on the sum of
                % edge weights to the previous region
                if r > 1
                    if sum(sum(G.W(rs0:rs1-1,rs1:rs1+n1-1))) < sum(sum(G.W(rs0:rs1-1,rs1+n1:rs2-1)))
                        pm = -pm;
                        n1 = n-n1;
                    end
                end
                
                % update the indexing
                ind(rs1:rs1+n1-1) = indrs(pm > 0);
                ind(rs1+n1:rs2-1) = indrs(pm < 0);

                % update the region tracking
                rs(rr+1,j+1) = rs1+n1;
                rs(rr+2,j+1) = rs2;
                rs0 = rs1+n1;
                rr = rr+2;
                
            % regions with 1 node
            elseif n == 1
                rs(rr+1,j+1) = rs2;
                rs0 = rs1;
                rr = rr+1;
            end
        end

        j = j+1;
    end
    
    
% L_rw
elseif nargin > 1
    while regioncount < N
        % the number of regions on level j
        regioncount = nnz(rs(:,j))-1;
        
        % add a column for level j+1, if necessary
        if j == jmax
            rs(1,j+1) = 1;
            jmax = jmax+1;
        end
        
        % for tracking the child regions
        rr = 1;
        
        % cycle through the parent regions
        for r = 1:regioncount
            % the start of the parent region and 1 node after the end of
            % the parent region
            rs1 = rs(r,j);
            rs2 = rs(r+1,j);
            
            % the number of nodes in the parent region
            n = rs2-rs1;
            
            % regions with 2 or more nodes
            if n > 1
                indrs = ind(rs1:rs2-1);

                % partition the current region
                pm = PartitionFiedlerSweep(G.W(indrs,indrs),1);
                
                % determine the number of points in child region 1
                n1 = sum(pm > 0);

                % switch regions 1 and 2, if necessary, based on the sum of
                % edge weights to the previous region
                if r > 1
                    if sum(sum(G.W(rs0:rs1-1,rs1:rs1+n1-1))) < sum(sum(G.W(rs0:rs1-1,rs1+n1:rs2-1)))
                        pm = -pm;
                        n1 = n-n1;
                    end
                end
                
                % update the indexing
                ind(rs1:rs1+n1-1) = indrs(pm > 0);
                ind(rs1+n1:rs2-1) = indrs(pm < 0);

                % update the region tracking
                rs(rr+1,j+1) = rs1+n1;
                rs(rr+2,j+1) = rs2;
                rs0 = rs1+n1;
                rr = rr+2;
                
            % regions with 1 node
            elseif n == 1
                rs(rr+1,j+1) = rs2;
                rs0 = rs1;
                rr = rr+1;
            end
        end

        j = j+1;
    end
end


% get rid of excess columns in rs
rs(:,j:end) = [];

% create a GraphPart object
if nargin == 1
    GP = GraphPart(ind,rs,[],[],[],[],[],'Fiedler (L)');
else
    GP = GraphPart(ind,rs,[],[],[],[],[],'Fiedler (L_rw)');
end


end