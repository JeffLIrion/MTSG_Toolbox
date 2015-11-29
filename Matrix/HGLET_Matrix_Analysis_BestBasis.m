function [dvec,BSrows,BScols] = HGLET_Matrix_Analysis_BestBasis(matrix,GProws,GPcols,Wrows,Wcols,costfun_rows,costfun_cols,flatten_rows,flatten_cols,~,~)
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
%   ~                   if 10 inputs are given, use the eigevectors of L_rw
%                       as the bases  ==>  U_rw = D^(0.5) * U_sym
%   ~                   if 11 inputs are given, use the eigenvectors of
%                       L_sym as the bases  ==>  U_sym = D^(-0.5) * U_rw
%
% Output
%   dvec                the HGLET expansion coefficients (not redundant)
%   BSrows              the best basis on the rows
%   BScols              the best basis on the columns
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
if nargin < 10
    dmatrix = HGLET_Analysis(Grows,GProws);
elseif nargin == 10
    dmatrix = HGLET_Analysis(Grows,GProws,1);
else
    dmatrix = HGLET_Analysis(Grows,GProws,1,1);
end

% analyze the data matrix using the columns
if nargin < 10
    dmatrix2 = HGLET_Analysis(Gcols,GPcols);
elseif nargin == 10
    dmatrix2 = HGLET_Analysis(Gcols,GPcols,1);
else
    dmatrix2 = HGLET_Analysis(Gcols,GPcols,1,1);
end

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
[dvec,BSrows] = HGLET_BestBasis(dmatrix,GProws,costfun_rows,flatten_rows);
[~,BScols] = HGLET_BestBasis(dmatrix2,GPcols,costfun_cols,flatten_cols);

% analyze the row best-basis coefficient matrix using the columns
if nargin < 10
    dmatrix = HGLET_Analysis(ReplaceData(Gcols,dvec'),GPcols);
elseif nargin == 10
    dmatrix = HGLET_Analysis(ReplaceData(Gcols,dvec'),GPcols,1);
else
    dmatrix = HGLET_Analysis(ReplaceData(Gcols,dvec'),GPcols,1,1);
end

% extract the col best-basis coefficients
dvec = dmatrix2dvec(dmatrix,GPcols,BScols)';


end