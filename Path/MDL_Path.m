function cost = MDL_Path(xQ,xQO,xO,q)
% Compute the quantitative MDL for a set of coefficients
%
% Input
%   xQ              the quantized coefficients <== divided by 2\delta, 
%                   round, and threshold
%   xQO             the quantized orthonormal coefficients ==> xQ but with
%                   Lrw coefficients modified for the purpose of computing
%                   2-norms
%   xO              the true orthonormal coefficients coefficients <== 
%                   divided by 2\delta, modify Lrw coefficients
%   q               the quantization precision is delta=2^-q
%   
% Output
%   cost            the quantitative MDL cost
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



n = length(xQ);

% assign Infinite cost to vectors of length 1
if n < 2
    cost = Inf;
    return
end

% define the precision (Inf-norm)
delta = 2^(-q);

% tabulate the integer coefficients and compute their distribution
K = max(abs(xQ));
p = tabulate([xQ;(-K:K)']);
p = p(:,2)-1;
p = p(p > 0)/n;

% the noise/residual
sigmasq = max(norm(2*delta*(xQO-xO),2)^2/n, 10^-13*0+eps);

% the quantization cost
cost = Huffman_cost(p,n) + Lstar(round(sigmasq/2/delta)) + n/2*log2(2*pi*exp(1)*sigmasq);


end




function y = Lstar(j)
% The codelength derived from the so-called "universal prior for integers"

c0 = 2.865064;

% j is a scalar
if length(j) == 1
    if j == 0
        y = 1;
    else
        y = log2(4*c0);
        logj = abs(j);
        while logj > 1
            logj = log2(logj);
            y = y+logj;
        end
    end
    
% j is a vector
else
    logj = abs(j);
    n = length(logj);
    nnzj = nnz(logj);
    
    % the zero entries in j
    y = n-nnzj;
    
    % the nonzero entries in j
    y = y+nnzj*log2(4*c0);
    logj = logj(logj > 1);
    while ~isempty(logj)
        logj = log2(logj);
        y = y+sum(logj);
        logj = logj(logj > 1);
    end
end
end




function y = Huffman_cost(p,N)
% The upper bound of the Huffman codelength of storing a vector of integers
% of length "N" with pmf "p" (where "p" contains no zero entries)

y = N*(1-dot(p,log2(p)));
end