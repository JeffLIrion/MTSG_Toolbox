function SNRs = nonorth2SNR(nonorth,f,GP,BS,SORH,B)
% Given a vector 'nonorth' of non-orthonormal expansion coefficients for a
% signal 'f', return a vector of SNR values when thresholding 1,2,...,N 
% coefficients.  
%
% Input
%   nonorth     a vector of non-orthonormal expansion coefficients
%   f           the noise-free signal
%   GP          a GraphPart object (used so the scaling coefficients aren't
%               thresholded)
%   BS          a BasisSpec object (used so the scaling coefficients aren't
%               thresholded)
%   SORH        use soft ('s') or hard ('h') thresholding
%   B           the matrix such that B*nonorth is the original noisy signal
%
% Output
%   SNRs        a vector of SNR values
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% constants
N = length(nonorth);
normf = norm(f,2);

% sort the coefficients
[~,IX] = sort(abs(nonorth),'descend');
nonorth = nonorth(IX);

% identify the scaling coefficients in the new sorted order
[~,~,tag] = ExtractData(GP);
tag = dmatrix2dvec(tag,GP,BS);
tag = double(tag == 0);
tag = tag(IX);    

% the coefficients to keep
Keep = double( triu(ones(N,N))+repmat(tag,1,N) > 0 );

% hard threshold
nonorthT = bsxfun(@times, nonorth, Keep);

% soft threshold
if strcmpi(SORH,'s') || strcmpi(SORH,'soft')
    ST = bsxfun(@times, [nonorth(2:end);0]', triu(repmat(double(tag==0),1,N)));
    nonorthT = nonorthT - sign(nonorthT).*abs(ST);
end

% put the coefficients back in their original order
nonorthT(IX,:) = nonorthT;

% compute the noise^2
noise = (sum(bsxfun(@minus, B*nonorthT, f).^2))';

% compute the SNRs
SNRs = 20*log10( normf * noise.^-0.5 );


end