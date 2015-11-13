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



% p-norms or p-quasinorms
if isnumeric(flatten)
    if flatten == 1
        dmatrix = sum(abs(dmatrix),3);
    elseif flatten == 0
        dmatrix = sum(dmatrix ~= 0,3);
    else
        dmatrix = (sum(abs(dmatrix).^flatten,3)).^(1/flatten);
    end
    
% histogram
% % % elseif strcmpi(flatten,'histogram')
    
    
% sum of absolute values = 1-norm
elseif strcmpi(flatten,'abs')
    dmatrix = sum(abs(dmatrix),3);
    
% standard deviation
elseif strcmpi(flatten,'std')
    dmatrix = std(dmatrix,0,3);
    
% variance
elseif strcmpi(flatten,'var')
    dmatrix = var(dmatrix,0,3);
    
% inverse variance
elseif strcmpi(flatten,'inverse variance') || strcmpi(flatten,'invvar')
    dmatrix = var(dmatrix,0,3);
    dmatrix( abs(dmatrix) < eps) = eps;
    dmatrix = dmatrix.^-1;
    
% max coefficient - min coefficient
elseif strcmpi(flatten,'minmax') || strcmpi(flatten,'maxmin')
    dmatrix = max(dmatrix,[],3)-min(dmatrix,[],3);
    
% sum (no absolute value)
elseif strcmpi(flatten,'sum')
    dmatrix = sum(dmatrix,3);
    
% the sum of the differences of consecutive coeffs
elseif strcmpi(flatten,'sumdiff')
    dmatrix = sum( dmatrix(:,:,2:end)-dmatrix(:,:,1:end-1), 3);
    
% the sum of the absolute values of the differences of consecutive coeffs
elseif strcmpi(flatten,'sumabsdiff')
    dmatrix = sum( abs( dmatrix(:,:,2:end)-dmatrix(:,:,1:end-1) ), 3);
    
% Shannon entropy
elseif strcmpi(flatten,'entropy')
    p = abs(dmatrix)/norm(dmatrix(:),'fro');
    dmatrix = sum(p.*log2(p),3);
    
% threshold coefficients below 0.5*norm(dmatrix)/numel(dmatrix) and sum
elseif strcmpi(flatten,'sub')
    t = 0.5*norm(dmatrix(:),'fro')/numel(dmatrix);
    dmatrix(abs(dmatrix) < t) = 0;
    dmatrix = sum(dmatrix,3);
    
% default (1-norm)
else
    dmatrix = sum(abs(dmatrix),3);
end


end