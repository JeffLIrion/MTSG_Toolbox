function [G,sigma] = AddNoise(G,snr,noisetype)
% Add noise to the data of a GraphSig object
%
% Input
%   G           a GraphSig object
%   snr         the SNR that the noisy signal should have
%   noisetype   the type of noise: Gaussian (default) or Poisson
%
% Outputs
%   G           the GraphSig object with added noise
%   sigma       the standard deviation of the noise
%
% 
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~exist('snr','var')
    snr = 1.24;
end

if ~exist('noisetype','var')
    noisetype = 'gaussian';
end


if strcmpi(noisetype,'gaussian')
    % generate Gaussian noise
    noise = randn(size(G.f));

    % scale the noise to the desired SNR level
    sigma = norm(G.f,'fro') / 10.0^(0.05*snr) / norm(noise,'fro');
    noise = sigma*noise;

    % generate the noisy signal
    G.f = G.f + noise;
    
elseif strcmpi(noisetype,'poisson')
    % generate Poisson noise
    f = double(full(G.f));
    noise = poissrnd(f)-f;
    
    % scale the noise to the desired SNR level
    sigma = norm(G.f,'fro') / 10.0^(0.05*snr) / norm(noise,'fro');
    noise = sigma*noise;
    
    % generate the noisy signal
    G.f = f + noise;
    
end

G.name = sprintf('%s (SNR = %.2f)',G.name,snr);

end