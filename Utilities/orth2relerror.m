function relerror = orth2relerror(orth)
% Given a vector 'orth' of orthonormal expansion coefficients, return a 
% vector of relative approximation errors when retaining the 1,2,...,N 
% largest coefficients in magnitude.  
%
% Input
%   orth        a vector of orthonormal expansion coefficients
%
% Output
%   relerror    a vector of relative approximation errors
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% sort the coefficients
orth = sort(orth.^2,'descend');

% compute the relative errors
relerror = ((abs(sum(orth)-cumsum(orth))).^0.5)/sum(orth).^0.5;


end