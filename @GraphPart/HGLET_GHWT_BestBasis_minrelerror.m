function [dvec,BS,trans,tau] = HGLET_GHWT_BestBasis_minrelerror(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GP,G,~)
% Find the best basis for approximating the signal 'G' by performing the
% best basis search with a range of tau-measures as cost functionals 
% (tau = 0.1, 0.2, ..., 1.9) and minimizing the relative error.  
%
% Input
%   dmatrixH    the matrix of HGLET expansion coefficients ==> eigenvectors
%               of L
%   dmatrixHrw  the matrix of HGLET expansion coefficients ==> eigenvectors
%               of Lrw
%   dmatrixHsym the matrix of HGLET expansion coefficients ==> eigenvectors
%               of Lsym
%   dmatrixG    the matrix of GHWT expansion coefficients
%   GP          a GraphPart object
%   G           the GraphSig object
%   ~           if a 7th input is given, don't compare the hybrid best
%               basis to the GHWT fine-to-coarse best basis
%
% Output
%   dvec        the vector of expansion coefficients corresponding to the 
%               bestbasis
%   BS          a BasisSpec object which specifies the best-basis
%   trans       specifies which transform was used for that portion of the
%               signal: 
%                   00 = HGLET with L
%                   01 = HGLET with Lrw
%                   10 = HGLET with Lsym
%                   11 = GHWT
%   tau         the tau that yields the smallest relative error
%
%
% 
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% cycle through various 'tau' values, find the best basis, and compare the
% sums of the relative errors
for tau_temp = 0.1:0.1:1.9
    % we are only considering the GHWT
    if ~isarray(dmatrixH) && ~isarray(dmatrixHrw) && ~isarray(dmatrixHsym) && isarray(dmatrixG)
        [dvec_temp,BS_temp] = GHWT_BestBasis(dmatrixG,GP,tau_temp);
        levlist = ExtractData(BS_temp);
        trans_temp = true(length(levlist),2);
        orthbasis = true;
        
    % we are considering 1 or more HGLET variations
    else
        [dvec_temp,BS_temp,trans_temp] = HGLET_GHWT_BestBasis(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GP,tau_temp);
    
        % check whether any HGLET Lrw basis vectors are in the best basis
        orthbasis = true;
        [rows,~] = size(trans_temp);
        for row = 1:rows
            if ~trans_temp(row,1) && trans_temp(row,2)
                orthbasis = false;
                break
            end
        end
    end
    
    % compute the relative errors
    if orthbasis
        relerror_temp = orth2relerror(dvec_temp);
    else
        B = HGLET_GHWT_Synthesis(eye(length(dvec_temp)),GP,BS_temp,trans_temp,G);
        relerror_temp = nonorth2relerror(dvec_temp,B);
    end
    sumrelerror_temp = sum(relerror_temp);
    
    % consider the GHWT fine-to-coarse best basis
    if nargin < 7 && isarray(dmatrixG)
        [~,~,~,~,dvec_f2c,BS_f2c] = GHWT_BestBasis(dmatrixG,GP,tau_temp);
        sumrelerror_f2c = sum(orth2relerror(dvec_f2c));
        if sumrelerror_f2c < sumrelerror_temp
            sumrelerror_temp = sumrelerror_f2c;
            dvec_temp = dvec_f2c;
            BS_temp = BS_f2c;
            levlist = ExtractData(BS_f2c);
            trans_temp = repmat([true,true],length(levlist),1);
        end
    end
    
    % compare to the current lowest sum of relative errors
    if ~exist('sumrelerror','var') || sumrelerror_temp < sumrelerror
        dvec = dvec_temp;
        BS = BS_temp;
        trans = trans_temp;
        sumrelerror = sumrelerror_temp;
        tau = tau_temp;
    end
end


end