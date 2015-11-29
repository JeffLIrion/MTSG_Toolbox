function [weightfun,use_sqrt,use_RegInvEuc] = weight_function_matrix(weightfun)
% Determine the weight function to be used for constructing hierarchical
% trees on the rows and columns of a matrix
%
% Inputs
%   weightfun   the specification for the weight function
%
% Outputs
%   weightfun   the weight function (as a function_handle)
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



use_sqrt = false;
use_RegInvEuc = false;

% weightfun is empty ==> use a Gaussian
if isempty(weightfun)
    weightfun = @(x,y) exp(-sum( bsxfun(@minus,x,y).^2, 2));
    
    
% weightfun is a function_handle ==> no action needed
elseif isa(weightfun,'function_handle')
    return
    
% weightfun is a string
elseif ischar(weightfun)
    % Gaussian
    if strcmpi(weightfun,'Gaussian')
        weightfun = @(x,y) exp(-sum( bsxfun(@minus,x,y).^2, 2));
        
    % apply a square root to l^2 distance
    elseif strcmpi(weightfun,'sqrt')
        weightfun = @(x,y) (sum( bsxfun(@minus,x,y).^2, 2)).^0.5;
        use_sqrt = true;
        
    % inverse Euclidean distance
    elseif strcmpi(weightfun,'Inverse Euclidean') || strcmpi(weightfun,'InvEuc')
        weightfun = @(x,y) (sum( bsxfun(@minus,x,y).^2, 2)).^(-0.5);
        
    % regularized inverse Euclidean distance
    elseif strcmpi(weightfun,'RegInvEuc')
        weightfun = @(x,y) (sum( bsxfun(@minus,x,y).^2, 2)).^(-0.5);
        use_RegInvEuc = true;
        
    % inverse earth mover's distance
    elseif strcmpi(weightfun,'InvEMD')
        weightfun = @(x,y) (sum( abs(bsxfun(@minus, cumsum(x), cumsum(y,2))), 2)).^(-1.0);
        
    % inverse normalized earth mover's distance
    elseif strcmpi(weightfun,'InvNormEMD')
        weightfun = @(x,y) (sum( abs(bsxfun(@minus, cumsumnorm(x), cumsumnorm(y))), 2)).^(-1.0);
        
    % string not recognized ==> use a Gaussian
    else
        weightfun = @(x,y) exp(-sum( bsxfun(@minus,x,y).^2, 2));
    end
    
% the class of weightfun is not recognized ==> use a Gaussian
else
    weightfun = @(x,y) exp(-sum( bsxfun(@minus,x,y).^2, 2));
end


end