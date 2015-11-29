function GP = PartialTreePathRatioCut(G,GP,BS,combine12)
% Generate a recursive partitioning for a path graph 'G' which includes the
% basis specified by 'BS' <==> 'GP'  
%
% Input
%   G           a GraphSig object
%   GP          a GraphPart object
%   BS          a BasisSpec object
%   combine12   'choose' means choose whether regions 1 and 2 get combined,
%               'yes' means regions 1 and 2 get combined,
%               'no' means regions 1 and 2 do not get combined
%
% Output
%   GP          a recursive partitioning which includes the input basis
%
%
% 
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% extract Data
W = ExtractData(G);
ind = ExtractData(GP);
[~,levlengths] = ExtractData(BS,GP);
levlengths = double(levlengths);

% constants
N = length(G);
regioncount = length(levlengths);
jmax = max(3*floor(log2(N)),4);

% stands for regionstarts, meaning that the index in 'ind' of the first
% point in region number i is rs(i)
sslac = ind_class(N);
rs = zeros(N+1,jmax,sslac);
rs(1,:) = 1;


%% 1) Proceed from the partitioning in BS to the coarsest level

% the regionstarts according to the partitioning in BS
rs(2:regioncount+1,1) = 1+cumsum(levlengths);

% combine pairs of regions
r = 1;
j = 1;
while regioncount > 1
    while r < regioncount
        rs2 = rs(r+1,j);
        rs3 = rs(r+2,j);
        
        % deal with regions 1 and 2
        if r == 1 && j == 1
            % 'choose' whether regions 1 and 2 should be combined
            if ~exist('combine12','var') || strcmpi(combine12,'choose')
                if regioncount == 2 || W(rs2-1,rs2) >= W(rs3-1,rs3)
                    levlengths(r) = levlengths(r)+levlengths(r+1);
                    levlengths(r+1) = [];
                    regioncount = regioncount-1;
                end
                
            % 'yes' ==> combine regions 1 and 2
            elseif strcmpi(combine12,'yes')
                levlengths(r) = levlengths(r)+levlengths(r+1);
                levlengths(r+1) = [];
                regioncount = regioncount-1;
                
            % 'no' ==> do not combine regions 1 and 2 (NO ACTION NEEDED)
            end
            
        % deal with regions other than 1 and 2
        else
            if r == regioncount-1 || W(rs2-1,rs2) >= W(rs3-1,rs3)
                levlengths(r) = levlengths(r)+levlengths(r+1);
                levlengths(r+1) = [];
                regioncount = regioncount-1;
            end
        end
        
        % move on to the next region
        r = r+1;
    end
    
    % prepare to go to the next coarser level
    regioncount = length(levlengths);
    rs(2:regioncount+1,j+1) = 1+cumsum(levlengths);
    r = 1;
    j = j+1;
end


%% 2) Proceed from the partitioning in BS to the finest level

% reverse the order so that the coarsest region is the first column
rs(:,1:j) = rs(:,j:-1:1);

% split regions into two via RatioCut
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
            indrs = ind(rs(r,j):rs(r+1,j)-1);

            % partition the current region
            pm = PartitionPathRatioCut(W(indrs,indrs));
            r1 = sum(pm > 0);% # points in subregion 1

            % update the region tracking
            rr = nnz(rs(:,j+1));
            rs(rr+1,j+1) = rs1+r1;
            rs(rr+2,j+1) = rs2+1;

        % regions with 1 node
        else
            rr = nnz(rs(:,j+1));
            rs(rr+1,j+1) = rs(r+1,j);
        end
    end

    % take it to the next level
    j = j+1;
end


% get rid of excess columns in rs
rs(:,j:end) = [];

GP = GraphPart(ind,rs);


end