function dmatrix = dmatrix_flatten(dmatrix,flatten)
% Flatten dmatrix using the method specified by the string "flatten"
%
% Input
%   dmatrix     the matrix of expansion coefficients
%   flatten     the method for flattening dmatrix
%
% Output
%   dmatrix     the flattened matrix
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if isnumeric(flatten)
    if flatten == 1
        dmatrix = sum(abs(dmatrix),3);
    elseif flatten == 0
        dmatrix = sum(dmatrix ~= 0,3);
    else
        dmatrix = (sum(abs(dmatrix).^flatten,3)).^(1/flatten);
    end
    
% elseif strcmpi(flatten,'histogram')
    
    
elseif strcmpi(flatten,'std')
    dmatrix = std(dmatrix,0,3);
    
elseif strcmpi(flatten,'var')
    dmatrix = var(dmatrix,0,3);
    
else
    dmatrix = sum(abs(dmatrix),3);
end


end