function sslac = tag_class(jmax)
% Given jmax, determine what class "tag" should be 
%
% Inputs
%   jmax        the number of levels in the recursive partitioning
%               (j = 1,...,jmax)
%
% Outputs
%   sslac       the class that tag should be
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if jmax-1 <= 8
    sslac = 'uint8';
elseif jmax-1 <= 16
    sslac = 'uint16';
elseif jmax-1 <= 32
    sslac = 'uint32';
else
    sslac = 'double';
end


end