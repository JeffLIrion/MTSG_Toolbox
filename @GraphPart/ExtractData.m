function [ind,rs,tag,compinfo,rsf2c,tagf2c,compinfof2c,method] = ExtractData(GP)
% Return the data in a GraphPart object
%
% Input
%   GP              the input GraphPart object whose data is to be 
%                   extracted
%
% Output
%   ind             ordering of the indices on the finest level
%   rs              regionstarts (coarse-to-fine) <==> the index in 'ind' 
%                   of the first point in region number i is rs(i)
%   tag             tag info for the GHWT coarse-to-fine basis
%   compinfo        indicates whether the coefficient was formed from 2
%                   coefficents (value is nonzero) or from only 1
%                   coefficient (value is zero); when a scaling and
%                   Haar-like coefficient are formed, their corresponding
%                   values in compinfo indicate the number of nodes in each
%                   of the 2 subregions
%   rsf2c           regionstarts (fine-to-coarse) <==> the index in 'ind' 
%                   of the first point in region number i is rs(i)
%   tagf2c          tag info for the GHWT coarse-to-fine basis
%   compinfof2c     indicates whether the coefficient was formed from 2
%                   coefficents (value is nonzero) or from only 1
%                   coefficient (value is zero); when a scaling and
%                   Haar-like coefficient are formed, their corresponding
%                   values in compinfo indicate the number of nodes in each
%                   of the 2 subregions
%   method          how the partition tree was constructed
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if nargout > 2 && isempty(GP.rsf2c)
    GP = GHWT_Info(GP);
end

ind = GP.ind;
rs = GP.rs;
tag = GP.tag;
compinfo = GP.compinfo;
rsf2c = GP.rsf2c;
tagf2c = GP.tagf2c;
compinfof2c = GP.compinfof2c;
method = GP.method;


end