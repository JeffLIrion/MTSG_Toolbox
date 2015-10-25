function costfun = cost_functional(costfun)
% Determine the cost functional to be used by the best-basis algorithm.  
%
% Inputs
%   costfun     the specification for the cost functional
%
% Outputs
%   costfun     the cost functional (as a function_handle)
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if isempty(costfun)
    costfun = @(x) norm(x,0.1);
elseif isnumeric(costfun)
    costfun = @(x) norm(x,costfun);
elseif ~isa(costfun,'function_handle')
    costfun = @(x) norm(x,0.1);
end


end