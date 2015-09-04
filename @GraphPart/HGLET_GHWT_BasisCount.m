function [count1,count2] = HGLET_GHWT_BasisCount(GP)
% Given a recursive graph partitioning, count the number of chooseable
% bases.
%
% Input
%   GP          a GraphPart object
%
% Output
%   count1      the total number of bases which can be selected from the
%               HGLET dictionary / GHWT coarse-to-fine dictionary
%   count2      the total number of bases which can be selected from the
%               GHWT fine-to-coarse dictionary, excluding those which are
%               in the coarse-to-fine dictionary
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% fill out GP
GP = GHWT_Info(GP);

% constants
[N,jmax] = size(GP.rs);
N = N-1;

% construct a matrix that will count the number of possible bases for each
% block
block_bases = ones(N,jmax);

% construct a logical version of the tag matrix, where indices of scaling
% coefficients are true and all other values are false
tag_logical = false(N+1,jmax);
tag_logical(:,jmax) = true;


% count the number of possible bases
for j = jmax-1:-1:1
    regioncount = nnz(GP.rs(:,j))-1;
    tag_logical(GP.rs(1:regioncount,j),j) = true;
    for r = 1:regioncount
        % the index that marks the start of the first subregion
        rs1 = GP.rs(r,j);

        % the index that is one after the end of the second subregion
        rs3 = GP.rs(r+1,j);

        % the number of points in the current region
        n = double(rs3 - rs1);

        if n > 1
            % the index that marks the start of the second subregion
            rs2 = rs1+1;
            while ~tag_logical(rs2,j+1) && rs2 < rs3 
                rs2 = rs2+1;
            end
            
            % the number of possible bases for the current region
            %%% case 1: it has a subregion
            if rs2 < rs3
                block_bases(rs1,j) = block_bases(rs1,j+1)*block_bases(rs2,j+1)+1;
                
            %%% case 2: it does not have a subregion
            else
                block_bases(rs1,j) = block_bases(rs1,j+1);
            end
        end
    end
end

count1 = block_bases(1,1);
count2 = count1-jmax;


end