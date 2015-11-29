function W = DualWeightMatrix(matrix,GP,weightfun)
% Compute the partition tree based weight matrix, or the "naive" pairwise
% weight matrix if no recursive partitioning is provided.  
%
% Input
%   matrix      the matrix whose columns are being used to generate W
%   GP          the partitioning tree on the rows of 'matrix'
%   weightfun   the general weight function
%
% Output
%   W           the partition-tree-based weight matrix if GP is a
%               proper GraphPart object OR the "naive" pairwise weight 
%               matrix if GP is not a GraphPart object
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)
% Wed Oct 14 21:06:07 2015  added by  Naoki Saito


    
% obtain constants
[N,fcols] = size(matrix);

% the weight function to be used
if ~exist('weightfun','var')
    weightfun = weight_function_matrix([]);
else
    weightfun = weight_function_matrix(weightfun);
end


%%% PARTITION TREE BASED WEIGHT MATRIX
if exist('GP','var') && isa(GP,'GraphPart')

    % extract data and obtain constants
    [ind,rs] = ExtractData(GP);
    [~,jmax] = size(rs);

    % the total number of regions in the recursive partitioning
    R = nnz(rs(2:end,:));

    % initialize d, where d is a matrix containing the averages over all
    % 'R' regions ==> each column of d contains a column from 'matrix' AND
    % the averages of that column over all (R-N) non-singleton regions
    d = [matrix; zeros(R-N,fcols)];
    
    % compute all the scaling coefficients
    dind = N+1;
    for j = jmax-1:-1:1
        regioncount = nnz(rs(:,j))-1;
        for r = 1:regioncount
            % the index that marks the start of the region
            rs1 = rs(r,j);

            % the index that is one after the end of the region
            rs2 = rs(r+1,j);

            % the number of points in the region
            n = double(rs2 - rs1);

            % the indices of the nodes in the region
            indrs = ind(rs1:rs2-1);

            % the scaling coefficient
            if n > 0
                d(dind,:) = sum(matrix(indrs,:),1)/n;
                dind = dind+1;
            else
                d(end,:) = [];
            end
        end
    end
    
    % construct the weight matrix
    W = squareform(pdist(d',weightfun));
    
    
%%% "NAIVE" PAIRWISE WEIGHT MATRIX
else
    W = squareform(pdist(matrix',weightfun));
end

% set the diagonal to 0
W(1:length(W)+1:numel(W)) = 0;


end