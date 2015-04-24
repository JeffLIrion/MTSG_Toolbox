function [dvec,BS,dmatrix,cmap] = HGLET_BestBasis(dmatrix,GP,costfun,flatten)
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


%% Find the best-basis

% allocate/initialize
dvec = dmatrix(:,jmax);
levlist = jmax*ones(N,1,'uint8');

% allocate a vector to store MDL costs
costs = zeros(N,1);
if useMDL
    for row = 1:N
        costs(row) = costfun(dvec(row));
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
            costBB = costfun(dvec(indr));
        end
        
        costNEW = costfun(dmatrix(indr,j));
        if costBB >= costNEW - tol
            [dvec(indr), levlist(indr), costs(indr)] = BBchange(costNEW,dmatrix(indr,j),j);
        end
    end
end

levlist( levlist==0 ) = [];

BS = BasisSpec(levlist,[],true,'HGLET Best Basis');
BS = levlist2levlengths(GP,BS);


% if using MDL, rescale the coefficients
if useMDL && fcols == 1
    dvec = dvec*dnorm;
end


if nargout > 2
    % define the color map -- color portion
    set(0, 'DefaultFigureVisible', 'off');
    cmap1 = colormap(jet);
    set(0, 'DefaultFigureVisible', 'on');
    
    % make the chosen basis in dmatrix negative so that it appears in 
    % grayscale and make the rest of the dictionary positive
    dmatrix = abs(dmatrix)+max(max(abs(dmatrix)))/length(cmap1);
    [~,levlengths] = ExtractData(BS);
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