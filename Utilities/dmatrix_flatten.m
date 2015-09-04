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
    
    
elseif strcmpi(flatten,'abs')
    dmatrix = sum(abs(dmatrix),3);
    
elseif strcmpi(flatten,'std')
    dmatrix = std(dmatrix,0,3);
    
elseif strcmpi(flatten,'var')
    dmatrix = var(dmatrix,0,3).^2;
    
elseif strcmpi(flatten,'minmax') || strcmpi(flatten,'maxmin')
    dmatrix = max(dmatrix,[],3)-min(dmatrix,[],3);
    
elseif strcmpi(flatten,'sum')
    dmatrix = sum(dmatrix,3);
    
elseif strcmpi(flatten,'sumdiff')
    dmatrix = sum( dmatrix(:,:,2:end)-dmatrix(:,:,1:end-1), 3);
    
elseif strcmpi(flatten,'sumabsdiff')
    dmatrix = sum( abs( dmatrix(:,:,2:end)-dmatrix(:,:,1:end-1) ), 3);
    
elseif strcmpi(flatten,'entropy')
    p = abs(dmatrix)/norm(dmatrix(:),'fro');
    dmatrix = sum(p.*log2(p),3);
    
elseif strcmpi(flatten,'sub')
    dmatrix = sqrt(numel(dmatrix))*abs(dmatrix)/norm(dmatrix(:),'fro');
    dmatrix(dmatrix < 0.5) = 0;
    dmatrix = sum(dmatrix,3);
    
else
    dmatrix = sum(abs(dmatrix),3);
end


end