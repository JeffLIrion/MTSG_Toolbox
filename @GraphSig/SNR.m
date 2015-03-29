function [value,sigma] = SNR(G1,G2)
% Compute the SNR between G1 (original signal) and G2 (noisy signal)
%
% Input
%   G1      Original reference signal
%   G2      Restored or noisy signal
%
% Output
%   value   Signal/Noise ratio
%   sigma   the standard deviation of the noise
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



value = 20*log10( norm(G1.f,'fro')/norm(G2.f-G1.f,'fro') );

if nargout == 2
    sigma = std(G2.f-G1.f);
end


end