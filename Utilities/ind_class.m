function sslac = ind_class(N)
% Given N, determine what class "ind" (and "rs") should be 
%
% Inputs
%   N           the length of the graph signal
%
% Outputs
%   sslac       the class that ind should be
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if N < 2^8-2
    sslac = 'uint8';
elseif N < 2^16-2
    sslac = 'uint16';
elseif N < 2^32-2
    sslac = 'uint32';
elseif N < 2^64-2
    sslac = 'uint64';
else
    sslac = 'double';
end


end