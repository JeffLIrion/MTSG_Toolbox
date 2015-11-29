function dvec = HGLET_Matrix_Analysis(matrix,GProws,GPcols,BSrows,BScols,Wrows,Wcols,~,~)
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
%   ~                   if 8 inputs are given, use the eigevectors of L_rw
%                       as the bases  ==>  U_rw = D^(0.5) * U_sym
%   ~                   if 9 inputs are given, use the eigenvectors of
%                       L_sym as the bases  ==>  U_sym = D^(-0.5) * U_rw
%
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
if nargin < 8
    dmatrix = HGLET_Analysis(Grows,GProws);
elseif nargin == 8
    dmatrix = HGLET_Analysis(Grows,GProws,1);
else
    dmatrix = HGLET_Analysis(Grows,GProws,1,1);
end

% extract the row best-basis coefficients
dvec = dmatrix2dvec(dmatrix,GProws,BSrows);

% analyze the row best-basis coefficient matrix using the columns
if nargin < 8
    dmatrix = HGLET_Analysis(ReplaceData(Gcols,dvec'),GPcols);
elseif nargin == 8
    dmatrix = HGLET_Analysis(ReplaceData(Gcols,dvec'),GPcols,1);
else
    dmatrix = HGLET_Analysis(ReplaceData(Gcols,dvec'),GPcols,1,1);
end

% extract the col best-basis coefficients
dvec = dmatrix2dvec(dmatrix,GPcols,BScols)';


end