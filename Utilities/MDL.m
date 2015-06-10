function cost = MDL(x,kmin,kmax,levlist_cost,levlens_cost,trans_cost)
% Compute the QMDL for a set of coefficients
%
% Input
%   x               the coefficient vector
%   kmin            the minimum k value <==> max delta = 2^-kmin
%   kmax            the maximum k value <==> min delta = 2^-kmax
%   levlist_cost    the cost of storing the level of the recursive
%                   partitioning from which a block of coefficients 
%                   originates
%   levlens_cost    the cost of storing the length of the current block of
%                   coefficients (unnecessary, since this is the same as
%                   the length of x?)
%   trans_cost      the cost of specifying the transform used for the
%                   current block of coefficients
%   
% Output
%   cost            the quantitative MDL cost
%
%
%
% Copyright 2014 Jeff Irion and Naoki Saito.  
%
% Implemented by Jeff Irion, 2015.  


n = length(x);

% compute the costs of delta + quantized coefficients + sigmasq "reward"
for k = kmin:kmax
    % define the precision
    delta = 2^(-k);

    % compute integer coefficients
    xint = round(x/2/delta);

    % compute sigmasq
    sigmasq = max(norm(x-2*delta*xint,2)^2/n, 10^-13);
    
    if k == kmin
        cost = k + Lstar(max(abs(xint))) + Huffman_cost(xint) + Lstar(round(sigmasq/2/delta)) + n/2*log2(2*pi*exp(1)*sigmasq);
    else
        costnew = k + Lstar(max(abs(xint))) + Huffman_cost(xint) + Lstar(round(sigmasq/2/delta)) + n/2*log2(2*pi*exp(1)*sigmasq);
        if costnew < cost
            cost = costnew;
        end
    end
end

cost = cost + levlist_cost + levlens_cost + trans_cost;
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




function y = Huffman_cost(j)
% The upper bound of the Huffman codelength of storing a vector j of
% integers

N = length(j);

% generate the pmf
p = tabulate([j;-1]);
p = p(2:end,2)/N;

% compute the Huffman codelength
y = N*(1-dot(p,log2(p)));
end