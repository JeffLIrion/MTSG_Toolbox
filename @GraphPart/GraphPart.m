function GP = GraphPart(ind,rs,tag,compinfo,rsf2c,tagf2c,compinfof2c,method)
% Constructor for a GraphPart object
%
% Input
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
% Output
%   GP              a struct with eight fields: ind, rs, tag, compinfo,
%                   rsf2c, tagf2c, compinfof2c, and method
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% ind
if exist('ind','var')
    datastruct.ind = ind;
else
    datastruct.ind = [];
end


%% rs
if exist('rs','var') && length(rs) == length(datastruct.ind)+1
    datastruct.rs = rs;
else
    datastruct.rs = [];
    datastruct.ind = [];
end


%% tag
if exist('tag','var') && norm(size(datastruct.rs)-size(tag)-[1,0]) == 0
    datastruct.tag = tag;
else
    datastruct.tag = [];
end


%% compinfo
if exist('compinfo','var') && norm(size(datastruct.rs)-size(compinfo)-[1,0]) == 0
    datastruct.compinfo = compinfo;
else
    datastruct.compinfo = [];
end


%% rsf2c
if exist('rsf2c','var') && norm(size(datastruct.rs)-size(rsf2c)) == 0
    datastruct.rsf2c = rsf2c;
else
    datastruct.rsf2c = [];
end


%% tagf2c
if exist('tagf2c','var') && norm(size(datastruct.tag)-size(tagf2c)) == 0
    datastruct.tagf2c = tagf2c;
else
    datastruct.tagf2c = [];
end


%% compinfof2c
if exist('compinfof2c','var') && norm(size(datastruct.compinfo)-size(compinfof2c)) == 0
    datastruct.compinfof2c = compinfof2c;
else
    datastruct.compinfof2c = [];
end


%% method
if exist('method','var') && ischar(method)
    datastruct.method = method;
else
    datastruct.method = 'unspecified';
end


GP = class(datastruct,'GraphPart');


end