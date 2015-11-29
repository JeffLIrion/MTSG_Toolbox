function GP = GP_Fix(GP,G)
% Given a "broken" recursive partitioning, in which regions are split into
% two or more subregions OR a region with multiple nodes is not 
% partitioned, "fix" it.  
%
% Input
%   GP      the "broken" recursive partitioning
%   G       the GraphSig object to which GP corresponds (provides W for us)
%
% Output
%   GP      the "fixed" recursive partitioning
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% extract data
[ind,rs] = ExtractData(GP);
W = ExtractData(G);

N = length(G);
[~,jmax] = size(rs);

% cycle through the levels
j = 1;
regioncount = 0;
while regioncount < N
    regioncount = nnz(rs(:,j))-1;
    
    % add a column for level j+1, if necessary
    if j == jmax
        rs(:,j+1) = rs(:,j);
        jmax = jmax+1;
    end
    
    for r = 1:regioncount
        % the index that marks the start of the region
        rs1 = rs(r,j);
        
        % the index that is one after the end of the region
        rs2 = rs(r+1,j);
        
        % the number of points in the region
        n = double(rs2 - rs1);
        
        % the index that marks the start of the first subregion
        IX_subrs1 = find(rs(:,j+1) == rs1,1,'last');
        
        % the index that is one after the end of the last subregion
        IX_subrs2 = find(rs(:,j+1) == rs2,1,'first');
        
        % the number of subregions
        nsubs = IX_subrs2 - IX_subrs1;
        
        %%% Problem 1: the region has 2 or more nodes but only 1 subregion
        %%%  Solution: look at finer and finer levels until there are 2 or
        %%%            more subregions
        if n > 1 && nsubs == 1
            jj = j;
            while nsubs == 1
                jj = jj+1;
                IX1 = find(rs(:,jj+1) == rs1,1,'last');
                IX2 = find(rs(:,jj+1) == rs2,1,'first');
                nsubs = IX2 - IX1;
            end
            rs(IX_subrs1+nsubs:N+1,j+1) = rs(IX_subrs2:N+2-nsubs,j+1);
            rs(IX_subrs1:IX_subrs1+nsubs,j+1) = rs(IX1:IX2,jj+1);
            IX_subrs2 = IX_subrs1+nsubs;
            for col = j+2:jj
                IX1b = find(rs(:,col) == rs1);
                IX2b = find(rs(:,col) == rs2);
                rs(IX1b+nsubs:N+1,col) = rs(IX2b:N+2-nsubs,col);
                rs(IX1b:IX1b+nsubs,col) = rs(IX1:IX2,jj+1);
            end
        end
        
        %%% Problem 2: the region has 3 or more subregions
        %%%  Solution: cluster the subregions into exactly 2 subregions
        if nsubs > 2
            for IX_cut = IX_subrs1:IX_subrs2-1
                indrs1 = ind(rs1:rs(IX_cut+1,j+1)-1);
                indrs2 = ind(rs(IX_cut+1,j+1):rs2-1);
                n1 = length(indrs1);
                n2 = length(indrs2);
                cutsum = sum(sum(W(indrs1,indrs2)))*(1/n1+1/n2);
                if IX_cut == IX_subrs1 || cutsum < mincut
                    mincut = cutsum;
                    IX_mincut = IX_cut;
                end
            end
            if j+1 == jmax
                rs(:,j+2) = rs(:,j+1);
                jmax = jmax+1;
            end
            rs(IX_subrs1+1,j+1) = rs(IX_mincut+1,j+1);
            rs(IX_subrs1+2:N+3-nsubs,j+1) = rs(IX_subrs2:N+1,j+1);
            rs(N+4-nsubs:N+1,j+1) = 0;
        end
    end
    j = j+1;
end

if sum(abs(rs(:,end)-rs(:,end-1))) == 0
    rs(:,end) = [];
end

GP = GraphPart(ind,rs);


end
