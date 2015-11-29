function relerror = nonorth2relerror(nonorth,B)
% Given a vector 'nonorth' of non-orthonormal expansion coefficients and 
% the matrix 'B' such that B*nonorth is the original signal, return a
% vector of relative approximation errors when retaining the 1,2,...,N 
% largest coefficients in magnitude.  
%
% Input
%   nonorth     a vector of non-orthonormal expansion coefficients
%   B           the matrix whose such that B*nonorth is the original signal
%
% Output
%   relerror    a vector of relative approximation errors
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% sort the expansion coefficients
[~,IX] = sort(abs(nonorth),'descend');

% generate a matrix where column j contains the j largest coefficients
matrix = triu(repmat(nonorth(IX),1,length(IX)));
matrix(IX,:) = matrix;

% the original signal and a matrix of reconstructions
f = B*nonorth;
recons = B*matrix;

% the relative errors
relerror = (sum( (bsxfun(@minus,f,recons)).^2,1 ).^0.5)'/norm(f,2);


end