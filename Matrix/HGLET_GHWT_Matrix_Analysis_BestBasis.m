function [dvec,BSrows,BScols,transrows,transcols] = HGLET_GHWT_Matrix_Analysis_BestBasis(matrix,GProws,GPcols,Wrows,Wcols,costfun_rows,costfun_cols,flatten_rows,flatten_cols)
% For a matrix, equipped with row and column weight recursive partitionings
% and weight matrices, generate the (non-redundant) matrix of HGLET
% expansion coefficients, using the best basis algorithm to select best
% bases for the rows and columns
%
% Input
%   matrix              the matrix to be analyzed
%   GProws              the recursive partitioning on the rows
%   GPcols              the recursive partitioning on the columns
%   Wrows               the weight matrix on the rows
%   Wcols               the weight matrix on the columns
%   costfun_rows        the cost functional to be used for the rows
%   costfun_cols        the cost functional to be used for the columns
%   flatten_rows        the method for flattening vector-valued data to
%                       scalar-valued data for the rows
%   flatten_cols        the method for flattening vector-valued data to
%                       scalar-valued data for the columns
%
% Output
%   dvec                the HGLET expansion coefficients (not redundant)
%   BSrows              the best basis on the rows
%   BScols              the best basis on the columns
%   transrows           specifies which transform was used for that portion
%                       of the signal (rows): 
%                           00 = HGLET with L
%                           01 = HGLET with Lrw
%                           10 = HGLET with Lsym
%                           11 = GHWT
%   transcols           specifies which transform was used for that portion
%                       of the signal (columns)
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% generate GraphSig objects for the matrix
Grows = GraphSig(Wrows,[],matrix);
Gcols = GraphSig(Wcols,[],matrix');

% analyze the data matrix using the rows
[dmatrixH,dmatrixHrw,dmatrixHsym] = HGLET_Analysis_All(Grows,GProws);
dmatrixG = GHWT_Analysis(Grows,GProws);

% analyze the data matrix using the columns
[dmatrixH2,dmatrixHrw2,dmatrixHsym2] = HGLET_Analysis_All(Gcols,GPcols);
dmatrixG2 = GHWT_Analysis(Gcols,GPcols);
    
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
[dvec,BSrows,transrows] = HGLET_GHWT_BestBasis(dmatrixH, dmatrixHrw, dmatrixHsym, dmatrixG, GProws, costfun_rows, flatten_rows);
[~,BScols,transcols]    = HGLET_GHWT_BestBasis(dmatrixH2,dmatrixHrw2,dmatrixHsym2,dmatrixG2,GPcols, costfun_cols, flatten_cols);

% analyze the row best-basis coefficient matrix using the columns
[dmatrixH,dmatrixHrw,dmatrixHsym] = HGLET_Analysis_All(ReplaceData(Gcols,dvec'),GPcols);
dmatrixG = GHWT_Analysis(ReplaceData(Gcols,dvec'),GPcols);

% extract the col best basis coefficients
dvec = dmatrices2dvec(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GPcols,BScols,transcols)';


end