function SNRs = orth2SNR(orth,orth_true,GP,BS,SORH)
% Given a vector 'orth' of orthonormal expansion coefficients for a signal 
% 'f', return a vector of SNR values  when thresholding 1,2,...,N 
% coefficients.  
%
% Input
%   orth        a vector of orthonormal expansion coefficients
%   orth_true   the noise-free orthonormal expansion coefficients
%   GP          a GraphPart object (used so the scaling coefficients aren't
%               thresholded)
%   BS          a BasisSpec object (used so the scaling coefficients aren't
%               thresholded)
%   SORH        use soft ('s') or hard ('h') thresholding
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
N = length(orth);
normf = norm(orth_true,2);

% sort the coefficients
[~,IX] = sort(abs(orth),'descend');
orth = orth(IX);
orth_true = orth_true(IX);

% identify the scaling coefficients in the new sorted order
[~,~,tag] = ExtractData(GP);
tag = dmatrix2dvec(tag,GP,BS);
tag = double(tag == 0);
tag = tag(IX);    

% the coefficients to keep
Keep = double( triu(ones(N,N))+repmat(tag,1,N) > 0 );

% hard threshold
orthT = bsxfun(@times, orth, Keep);

% soft threshold
if strcmpi(SORH,'s') || strcmpi(SORH,'soft')
    ST = bsxfun(@times, [orth(2:end);0]', triu(repmat(double(tag==0),1,N)));
    orthT = orthT - sign(orthT).*abs(ST);
end

% compute the noise^2
noise = (sum( bsxfun(@minus, orthT, orth_true).^2 ))';

% compute the SNRs
SNRs = 20*log10( normf * noise.^-0.5 );


end