function [dvec,BS,dvecc2f,BSc2f,dvecf2c,BSf2c,dmatrix,cmap] = GHWT_BestBasis(dmatrix,GP,costfun,flatten)
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
%   dmatrix     a modified version of the coefficient matrix, with non-best
%               entries being positive and best basis entries being
%               negative
%   cmap        the color map for plotting dmatrix
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
[costfun,useMDL] = cost_functional(costfun);

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

% MDL stuff
if useMDL
    % compute the number of bits needed to store trans entries
    trans_cost = 0;
    
    % compute the number of bits needed to store levlist entries
    levlist_cost = 0*ceil(log2(jmax));
    
    % compute the number of bits needed to store levlengths entries
    % (equivalently, to store regionstarts entries)
    levlens_cost = ceil(log2(N));
    
    kmin = 1;
    kmax = 1 + ceil(0.5*log2(N));
    
    % define the cost functional
    costfun = @(x) MDL(x,kmin,kmax,levlist_cost,levlens_cost,trans_cost);
    
    % normalize the coefficients
    dnorm = norm(dmatrix(:,end),2)/sqrt(N);
    dmatrix = dmatrix/dnorm;
end


%% Find the best-basis from the GHWT coarse-to-fine dictionary

% allocate/initialize
dvecc2f = dmatrix(:,jmax);
levlistc2f = jmax*ones(N,1,'uint8');

% allocate a vector to store MDL costs
costs = zeros(N,1);
if useMDL
    for row = 1:N
        costs(row) = costfun(dvecc2f(row));
    end
end

% set the tolerance
tol = 10^4*eps;

% perform the basis search
for j = jmax-1:-1:1
    regioncount = nnz(GP.rs(:,j))-1;
    for r = 1:regioncount
        indr = GP.rs(r,j):GP.rs(r+1,j)-1;
        %%%%% compute the cost of the current best basis
        if useMDL
            costBB = sum(costs(indr));
        else
            costBB = costfun(dvecc2f(indr));
        end
        
        costNEW = costfun(dmatrix(indr,j));
        if costBB >= costNEW - tol
            [dvecc2f(indr), levlistc2f(indr), costs(indr)] = BBchange(costNEW,dmatrix(indr,j),j);
        end
    end
end

levlistc2f( levlistc2f==0 ) = [];

BSc2f = BasisSpec(levlistc2f,[],true,'GHWT c2f Best Basis');
BSc2f = levlist2levlengths(GP,BSc2f);
costc2f = costfun(dvecc2f);


% if using MDL, rescale the coefficients
if useMDL && fcols == 1
    dvecc2f = dvecc2f*dnorm;
    costc2f = sum(costs);
end


%% Find the best-basis from the GHWT fine-to-coarse dictionary

% generate the fine-to-coarse GraphPart fields and coefficient matrix
[GP,dmatrixf2c] = FineToCoarse(GP,dmatrix);

% allocate/initialize
dvecf2c = dmatrixf2c(:,jmax);
levlistf2c = jmax*ones(N,1,'uint8');

% allocate a vector to store MDL costs
costs = zeros(N,1);
if useMDL
    for row = 1:N
        costs(row) = costfun(dvecf2c(row));
    end
end

% perform the basis search
for j = jmax-1:-1:1
    regioncount = nnz(GP.rsf2c(:,j))-1;
    for r = 1:regioncount
        indr = GP.rsf2c(r,j):GP.rsf2c(r+1,j)-1;
        %%%%% compute the cost of the current best basis
        if useMDL
            costBB = sum(costs(indr));
        else
            costBB = costfun(dvecf2c(indr));
        end
        
        costNEW = costfun(dmatrixf2c(indr,j));
        if costBB >= costNEW - tol
            [dvecf2c(indr), levlistf2c(indr), costs(indr)] = BBchange(costNEW,dmatrixf2c(indr,j),j);
        end
    end
end

levlistf2c( levlistf2c==0 ) = [];

BSf2c = BasisSpec(levlistf2c,[],false,'GHWT f2c Best Basis');
BSf2c = levlist2levlengths(GP,BSf2c);
costf2c = costfun(dvecf2c);


% if using MDL, rescale the coefficients
if useMDL && fcols == 1
    dvecf2c = dvecf2c*dnorm;
    costf2c = sum(costs);
end


%% Compare the coarse-to-fine and fine-to-coarse best-bases

if costc2f >= costf2c
    dvec = dvecf2c;
    BS = BSf2c;
    if nargout > 6
        dmatrix = dmatrixf2c;
    end
else
    dvec = dvecc2f;
    BS = BSc2f;
end


if nargout > 6
    % define the color map -- color portion
    set(0, 'DefaultFigureVisible', 'off');
    cmap1 = colormap(jet);
    set(0, 'DefaultFigureVisible', 'on');
    
    % make the chosen basis in dmatrix negative so that it appears in 
    % grayscale and make the rest of the dictionary positive
    dmatrix = abs(dmatrix)+max(max(abs(dmatrix)))/length(cmap1);
    [levlist,levlengths] = ExtractData(BS);
    n = 1;
    for row = 1:length(levlist)
        dmatrix(n:n+levlengths(row)-1,levlist(row),:) = -abs(dvec(n:n+levlengths(row)-1,:));
        n = n+levlengths(row);
    end
    
    % define the color map -- grayscale portion
    dmax = max(max(dmatrix));
    dmin = max(abs(dvec));
    numgray = ceil(dmin*length(cmap1)/dmax);
    cmap2 = (linspace(0,1,numgray))';
    cmap = [cmap2, cmap2, cmap2; cmap1];
end


% if we flattened dmatrix, then "unflatten" the expansion coefficients
if fcols > 1
    dvec = dmatrix2dvec(dmatrix0,GP,BS);
end


end




function [dvec, levlist, costs] = BBchange(costNEW, dvec,j)
% Change to the new best basis

n = length(dvec);

levlist = zeros(n,1,'uint8');
levlist(1) = j;

costs = zeros(n,1);
costs(1) = costNEW;
end