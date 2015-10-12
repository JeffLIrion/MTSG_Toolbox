function TF = isarray(x)
% Return true if x is a non-empty numerical array and false otherwise
% 
% Input
%   x       a variable
% 
% Output
%   TF      true if x is a non-empty numerical array and false otherwise
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~isempty(x) && isnumeric(x) && ~isscalar(x)
    TF = true;
else
    TF = false;
end


end