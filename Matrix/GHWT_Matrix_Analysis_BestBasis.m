function [dvec,BSrows,BScols] = GHWT_Matrix_Analysis_BestBasis(matrix,GProws,GPcols,costfun_rows,costfun_cols,flatten_rows,flatten_cols)
% For a matrix, equipped with row and column weight recursive partitionings
% and weight matrices, generate the (non-redundant) matrix of GHWT
% expansion coefficients, using the best basis algorithm to select best
% bases for the rows and columns
%
% Input
%   matrix              the matrix to be analyzed
%   GProws              the recursive partitioning on the rows
%   GPcols              the recursive partitioning on the columns
%   costfun_rows        the cost functional to be used for the rows
%   costfun_cols        the cost functional to be used for the columns
%   flatten_rows        the method for flattening vector-valued data to
%                       scalar-valued data for the rows
%   flatten_cols        the method for flattening vector-valued data to
%                       scalar-valued data for the columns
% 
% Output
%   dvec                the GHWT expansion coefficients (not redundant)
%   BSrows              the best basis on the rows
%   BScols              the best basis on the columns
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

% analyze the data matrix using the columns
dmatrix2 = GHWT_Analysis(Gcols,GPcols);

% use default costfun and flatten values if they aren't provided
if ~exist('costfun_rows','var')
    costfun_rows = 1;
end
if ~exist('costfun_cols','var')
    costfun_cols = 1;
end
if ~exist('flatten_rows','var')
    flatten_rows = 1;
end
if ~exist('flatten_cols','var')
    flatten_cols = 1;
end

% find the row and column best bases
[dvec,BSrows] = GHWT_BestBasis(dmatrix, GProws,costfun_rows,flatten_rows);
[~,BScols] = GHWT_BestBasis(dmatrix2,GPcols,costfun_cols,flatten_cols);

% analyze the row best-basis coefficient matrix using the columns
dmatrix = GHWT_Analysis(ReplaceData(Gcols,dvec'),GPcols);

% extract the col best-basis coefficients
dvec = dmatrix2dvec(dmatrix,GPcols,BScols)';


end