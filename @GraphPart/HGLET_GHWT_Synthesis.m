function [f,GS] = HGLET_GHWT_Synthesis(dvec,GP,BS,trans,G)
% Given a vector of HGLET & GHWT expansion coefficients, info about the 
% graph partitioning, and the choice of basis and corresponding transforms,
% reconstruct the signal
%
% Input
%   dvec        the expansion coefficients corresponding to the chosen
%               basis
%   GP          a GraphPart object
%   BS          a BasisSpec object
%   trans       specifies which transform was used for that portion of the
%               signal: 
%                   00 = HGLET with L
%                   01 = HGLET with Lrw
%                   10 = HGLET with Lsym
%                   11 = GHWT
%   G           a GraphSig object
%
% Output
%   f           the reconstructed signal
%   GS          the reconstructed GraphSig object
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% fill out trans
[~,~,transfull] = BSfull(GP,BS,trans);

% decompose dvec into GHWT and HGLET components
dvecH    = bsxfun(@times, dvec, ~transfull(:,1) .* ~transfull(:,2));
dvecHrw  = bsxfun(@times, dvec, ~transfull(:,1) .*  transfull(:,2));
dvecHsym = bsxfun(@times, dvec,  transfull(:,1) .* ~transfull(:,2));
dvecG    = bsxfun(@times, dvec,  transfull(:,1) .*  transfull(:,2));


%% Synthesize using the transforms separately

fH = HGLET_Synthesis(dvecH,GP,BS,G);
fHrw = HGLET_Synthesis(dvecHrw,GP,BS,G,1);
fHsym = HGLET_Synthesis(dvecHsym,GP,BS,G,1,1);
fG = GHWT_Synthesis(dvecG,GP,BS);

f = fH + fHrw + fHsym + fG;

% create a GraphSig object with the reconstructed data
if nargout == 2
    if exist('G','var') == 1
        GS = ReplaceData(G,f);
    else
        GS = 0;
    end
end


end