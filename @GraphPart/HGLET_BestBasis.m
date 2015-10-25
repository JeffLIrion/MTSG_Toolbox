function [dvec,BS] = HGLET_BestBasis(dmatrix,GP,costfun,flatten)
% Select the best basis from the matrix of HGLET expansion coefficients
%
% Input
%   dmatrix     the matrix of expansion coefficients
%   GP          a GraphPart object
%   costfun     the cost functional to be used
%   flatten     the method for flattening vector-valued data to
%               scalar-valued data
%
% Output
%   dvec        the vector of expansion coefficients corresponding to the 
%               bestbasis
%   BS          a BasisSpec object which specifies the best-basis
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% the cost functional to be used
if ~exist('costfun','var')
    costfun = [];
end
costfun = cost_functional(costfun);

% constants and dmatrix cleanup
[N,jmax,fcols] = size(dmatrix);
dmatrix( abs(dmatrix) < 10^2*eps ) = 0;

% "flatten" dmatrix
if fcols > 1
    dmatrix0 = dmatrix;
    if ~exist('flatten','var')
        flatten = 1;
    end
    dmatrix = dmatrix_flatten(dmatrix,flatten);
end


%% Find the best-basis

% allocate/initialize
dvec = dmatrix(:,jmax);
levlist = jmax*ones(N,1,'uint8');

% set the tolerance
tol = 10^4*eps;

% perform the basis search
for j = jmax-1:-1:1
    regioncount = nnz(GP.rs(:,j))-1;
    for r = 1:regioncount
        indr = GP.rs(r,j):GP.rs(r+1,j)-1;
        %%%%% compute the cost of the current best basis
        costBB = costfun(dvec(indr));
        
        costNEW = costfun(dmatrix(indr,j));
        if costBB >= costNEW - tol
            [dvec(indr), levlist(indr)] = BBchange(dmatrix(indr,j),j);
        end
    end
end

levlist( levlist==0 ) = [];

BS = BasisSpec(levlist,[],true,'HGLET Best Basis');
BS = levlist2levlengths(GP,BS);


% if we flattened dmatrix, then "unflatten" the expansion coefficients
if fcols > 1
    dvec = dmatrix2dvec(dmatrix0,GP,BS);
end


end




function [dvec, levlist] = BBchange(dvec,j)
% Change to the new best basis

n = length(dvec);

levlist = zeros(n,1,'uint8');
levlist(1) = j;
end