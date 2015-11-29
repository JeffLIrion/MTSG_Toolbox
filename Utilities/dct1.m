function d = dct1(f)
% Perform the DCT-I transform on the column(s) of f.  
%
% Input
%   f       a column vector or matrix of column vectors
%
% Output
%   d       the DCT-I transform of (the columns of) f
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if length(f) == 1
    d = f;
    return
end

[N,~] = size(f);
N = N-1;
FUf = ifft([sqrt(2)*f(1,:); f(2:end-1,:); sqrt(2)*f(end,:); f(end-1:-1:2,:)]);
d = [sqrt(2)*FUf(1,:); FUf(2:N,:) + FUf(end:-1:N+2,:); sqrt(2)*FUf(N+1,:)]*sqrt(N/2);


end