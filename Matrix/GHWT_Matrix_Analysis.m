function dvec = GHWT_Matrix_Analysis(matrix,GProws,GPcols,BSrows,BScols)
% For a matrix, equipped with row and column weight recursive partitionings
% and weight matrices, generate the (non-redundant) matrix of GHWT
% expansion coefficients corresponding to the tensor product of specified
% row and column bases
%
% Input
%   matrix              the matrix to be analyzed
%   GProws              the recursive partitioning on the rows
%   GPcols              the recursive partitioning on the columns
%   BSrows              the best basis on the rows
%   BScols              the best basis on the columns
% 
% Output
%   dvec                the GHWT expansion coefficients (not redundant)
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% generate GraphSig objects for the matrix
[rows,cols] = size(matrix);
Grows = GraphSig(sparse(rows,rows),[],matrix);
Gcols = GraphSig(sparse(cols,cols),[],matrix');

% analyze the data matrix using the rows
dmatrix = GHWT_Analysis(Grows,GProws);

% extract the row best-basis coefficients
dvec = dmatrix2dvec(dmatrix,GProws,BSrows);

% analyze the row best-basis coefficient matrix using the columns
dmatrix = GHWT_Analysis(ReplaceData(Gcols,dvec'),GPcols);

% extract the col best-basis coefficients
dvec = dmatrix2dvec(dmatrix,GPcols,BScols)';


end