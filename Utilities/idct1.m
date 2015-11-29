function f = idct1(d)
% Perform the inverse DCT-I transform on the column(s) of d.  
%
% Input
%   d       DCT-I expansion coefficients (either a column vector or a 
%           matrix of column vectors)
%
% Output
%   f       the reconstructed signal
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if length(d) == 1
    f = d;
    return
end

[N,~] = size(d);
N = N-1;
FUd = fft([sqrt(2)*d(1,:); d(2:end-1,:); sqrt(2)*d(end,:); d(end-1:-1:2,:)]);
f = [sqrt(2)*FUd(1,:); FUd(2:N,:) + FUd(end:-1:N+2,:); sqrt(2)*FUd(N+1,:)]/sqrt(8*N);


end