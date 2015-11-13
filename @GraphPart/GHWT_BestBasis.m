function [dvec,BS,dvecc2f,BSc2f,dvecf2c,BSf2c] = GHWT_BestBasis(dmatrix,GP,costfun,flatten)
% Select the best basis from the matrix of GHWT expansion coefficients
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
%               best basis
%   BS          a BasisSpec object which specifies the best basis
%   dvecc2f     the vector of expansion coefficients corresponding to the 
%               coarse-to-fine best basis
%   BSc2f       a BasisSpec object which specifies the coarse-to-fine best
%               basis
%   dvecf2c     the vector of expansion coefficients corresponding to the 
%               fine-to-coarse best basis
%   BSf2c       a BasisSpec object which specifies the fine-to-coarse best
%               basis
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


%% Find the best-basis from the GHWT coarse-to-fine dictionary

% allocate/initialize
dvecc2f = dmatrix(:,jmax);
levlistc2f = jmax*ones(N,1,'uint8');

% set the tolerance
tol = 10^4*eps;

% perform the basis search
for j = jmax-1:-1:1
    regioncount = nnz(GP.rs(:,j))-1;
    for r = 1:regioncount
        indr = GP.rs(r,j):GP.rs(r+1,j)-1;
        %%%%% compute the cost of the current best basis
        costBB = costfun(dvecc2f(indr));
        
        costNEW = costfun(dmatrix(indr,j));
        if costBB >= costNEW - tol
            [dvecc2f(indr), levlistc2f(indr)] = BBchange(dmatrix(indr,j),j);
        end
    end
end

levlistc2f( levlistc2f==0 ) = [];

BSc2f = BasisSpec(levlistc2f,[],true,'GHWT c2f Best Basis');
BSc2f = levlist2levlengths(GP,BSc2f);
costc2f = costfun(dvecc2f);


%% Find the best-basis from the GHWT fine-to-coarse dictionary

% generate the fine-to-coarse GraphPart fields and coefficient matrix
[GP,dmatrixf2c] = FineToCoarse(GP,dmatrix);

% allocate/initialize
dvecf2c = dmatrixf2c(:,jmax);
levlistf2c = jmax*ones(N,1,'uint8');

% perform the basis search
for j = jmax-1:-1:1
    regioncount = nnz(GP.rsf2c(:,j))-1;
    for r = 1:regioncount
        indr = GP.rsf2c(r,j):GP.rsf2c(r+1,j)-1;
        %%%%% compute the cost of the current best basis
        costBB = costfun(dvecf2c(indr));
        
        costNEW = costfun(dmatrixf2c(indr,j));
        if costBB >= costNEW - tol
            [dvecf2c(indr), levlistf2c(indr)] = BBchange(dmatrixf2c(indr,j),j);
        end
    end
end

levlistf2c( levlistf2c==0 ) = [];

BSf2c = BasisSpec(levlistf2c,[],false,'GHWT f2c Best Basis');
BSf2c = levlist2levlengths(GP,BSf2c);
costf2c = costfun(dvecf2c);


%% Compare the coarse-to-fine and fine-to-coarse best-bases

if costc2f >= costf2c
    dvec = dvecf2c;
    BS = BSf2c;
else
    dvec = dvecc2f;
    BS = BSc2f;
end


% if we flattened dmatrix, then "unflatten" the expansion coefficients
if fcols > 1
    dvec = dmatrix2dvec(dmatrix0,GP,BS);
    if nargout > 2
        dvecc2f = dmatrix2dvec(dmatrix0,GP,BSc2f);
        dvecf2c = dmatrix2dvec(dmatrix0,GP,BSf2c);
    end
end


end




function [dvec, levlist] = BBchange(dvec,j)
% Change to the new best basis

n = length(dvec);

levlist = zeros(n,1,'uint8');
levlist(1) = j;
end