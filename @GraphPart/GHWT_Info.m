function GP = GHWT_Info(GP)
% For a GraphPart object GP, use rs and ind to compute tag, compinfo,
% rsf2c, tagf2c, and compinfof2c
%
% Input
%   GP      the GraphPart object, with rs and ind
%
% Output
%   GP      the GraphPart object, with all fields filled in
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% fill in "tag" and "compinfo"
if isempty(GP.tag) || isempty(GP.compinfo)
    GP = GHWT_Core(GP);
end

% fill in the fine-to-coarse fields
GP = FineToCoarse(GP);


end