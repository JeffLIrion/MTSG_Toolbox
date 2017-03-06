function BS = LevelBasisSpec(GP,j,~)
% Specify the basis corresponding to level j for a given graph partitioning
%
% Input
%   GP      a GraphPart object
%   j       the level to which the basis corresponds (j=0 is the global
%           level)
%   ~       if a 3rd input is given, use the fine-to-coarse dictionary
% 
% Output
%   BS      a BasisSpec object corresponding to the level j basis
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% coarse-to-fine dictionary
if nargin == 2
    Nj = nnz(GP.rs(:,j+1))-1;
    levlist = (j+1)*ones(Nj,1,'uint8');
    BS = BasisSpec(levlist,[],true,sprintf('coarse-to-fine level %d',j));
    
% fine-to-coarse dictionary
elseif nargin == 3
    if isempty(GP.rsf2c)
        GP = GHWT_Info(GP);
    end
    Nj = nnz(GP.rsf2c(:,j+1))-1;
    levlist = (j+1)*ones(Nj,1,'uint8');
    BS = BasisSpec(levlist,[],false,sprintf('fine-to-coarse level %d',j));
end

BS = levlist2levlengths(GP,BS);


end
