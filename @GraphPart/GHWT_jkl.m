function [j,k,l] = GHWT_jkl(GP,drow,dcol)
% Generate the (j,k,l) indices for the GHWT basis vector corresponding to 
% the coefficient dmatrix(drow,dcol).  
%
% Inputs
%   GP          a GraphPart object
%   drow        the row of the expansion coefficient
%   dcol        the column of the expansion coefficient
%
% Outputs
%   j           the level index of the expansion coefficient
%   k           the subregion index of the expansion coefficient
%   l           the tag of the expansion coefficient
% 
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



j = dcol-1;

k = find( GP.rs(:,dcol) > drow, 1);
k = k-2;

if isempty(GP.tag)
    GP = GHWT_Core(GP);
end

l = GP.tag(drow,dcol);


end