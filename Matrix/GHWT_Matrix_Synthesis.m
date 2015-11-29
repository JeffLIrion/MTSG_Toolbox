function matrix = GHWT_Matrix_Synthesis(dvec,GProws,GPcols,BSrows,BScols)
% Given a vector of HGLET expansion coefficients and info about the graph 
% partitioning and the choice of basis, reconstruct the signal
%
% Input
%   dvec        the expansion coefficients corresponding to the chosen
%               basis
%   GProws      the recursive partitioning on the rows
%   GPcols      the recursive partitioning on the columns
%   BSrows      the best basis on the rows
%   BScols      the best basis on the columns
%
% Output
%   matrix      the reconstructed matrix
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% reconstruct
dvec = GHWT_Synthesis(dvec',GPcols,BScols);
matrix = GHWT_Synthesis(dvec',GProws,BSrows);


end