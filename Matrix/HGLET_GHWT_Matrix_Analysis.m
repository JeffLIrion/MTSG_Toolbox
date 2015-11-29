function dvec = HGLET_GHWT_Matrix_Analysis(matrix,GProws,GPcols,BSrows,BScols,Wrows,Wcols,transrows,transcols)
% For a matrix, equipped with row and column weight recursive partitionings
% and weight matrices, generate the (non-redundant) matrix of HGLET
% expansion coefficients corresponding to the tensor product of specified
% row and column bases
%
% Input
%   matrix              the matrix to be analyzed
%   GProws              the recursive partitioning on the rows
%   GPcols              the recursive partitioning on the columns
%   BSrows              the best basis on the rows
%   BScols              the best basis on the columns
%   Wrows               the weight matrix on the rows
%   Wcols               the weight matrix on the columns
%   transrows           specifies which transform was used for that portion
%                       of the signal (rows): 
%                           00 = HGLET with L
%                           01 = HGLET with Lrw
%                           10 = HGLET with Lsym
%                           11 = GHWT
%   transcols           specifies which transform was used for that portion
%                       of the signal (columns)
% Output
%   dvec                the HGLET expansion coefficients (not redundant)
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

% extract the row best basis coefficients
dvec = dmatrices2dvec(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GProws,BSrows,transrows);

% analyze the row best-basis coefficient matrix using the columns
[dmatrixH,dmatrixHrw,dmatrixHsym] = HGLET_Analysis_All(ReplaceData(Gcols,dvec'),GPcols);
dmatrixG = GHWT_Analysis(ReplaceData(Gcols,dvec'),GPcols);

% extract the col best basis coefficients
dvec = dmatrices2dvec(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GPcols,BScols,transcols)';


end