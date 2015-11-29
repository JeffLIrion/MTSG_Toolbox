function dvec = dmatrices2dvec(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GP,BS,trans)
% Given matrices of HGLET & GHWT expansion coefficients, info about the 
% graph partitioning, and specifications for the basis and transforms used,
% return a vector of hybrid coefficients.  
%
% Input
%   dmatrixH        the matrix of HGLET expansion coefficients for L
%   dmatrixHrw      the matrix of HGLET expansion coefficients for Lrw
%   dmatrixHsym     the matrix of HGLET expansion coefficients for Lsym
%   dmatrixG        the matrix of GHWT expansion coefficients
%   GP              a GraphPart object
%   BS              a BasisSpec object
%   trans           specifies which transform was used for that portion of
%                   the signal: 
%                       00 = HGLET with L
%                       01 = HGLET with Lrw
%                       10 = HGLET with Lsym
%                       11 = GHWT
%
% Output
%   dvec            the vector of expansion coefficients corresponding to 
%                   the specified basis
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% fill out trans
[~,~,transfull] = BSfull(GP,BS,trans);

if isarray(dmatrixH)
    dvecH    = dmatrix2dvec(dmatrixH,GP,BS);
    dvecH    = bsxfun(@times, dvecH, ~transfull(:,1) .* ~transfull(:,2));
else
    dvecH    = 0;
end

if isarray(dmatrixHrw)
    dvecHrw  = dmatrix2dvec(dmatrixHrw,GP,BS);
    dvecHrw  = bsxfun(@times, dvecHrw, ~transfull(:,1) .*  transfull(:,2));
else
    dvecHrw  = 0;
end

if isarray(dmatrixHsym)
    dvecHsym = dmatrix2dvec(dmatrixHsym,GP,BS);
    dvecHsym = bsxfun(@times, dvecHsym, transfull(:,1) .* ~transfull(:,2));
else
    dvecHsym = 0;
end

if isarray(dmatrixG)
    dvecG    = dmatrix2dvec(dmatrixG,GP,BS);
    dvecG    = bsxfun(@times, dvecG, transfull(:,1) .*  transfull(:,2));
else
    dvecG    = 0;
end

dvec = dvecH + dvecHrw + dvecHsym + dvecG;


end