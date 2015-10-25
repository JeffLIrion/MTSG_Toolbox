function [dvec,BS,trans] = HGLET_GHWT_BestBasis(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GP,costfun,flatten)
% Select the best basis from several matrices of expansion coefficients
%
% Input
%   dmatrixH    the matrix of HGLET expansion coefficients ==> eigenvectors
%               of L
%   dmatrixHrw  the matrix of HGLET expansion coefficients ==> eigenvectors
%               of Lrw
%   dmatrixHsym the matrix of HGLET expansion coefficients ==> eigenvectors
%               of Lsym
%   dmatrixG    the matrix of GHWT expansion coefficients
%   GP          a GraphPart object
%   costfun     the cost functional to be used
%   flatten     the method for flattening vector-valued data to
%               scalar-valued data
%
% Output
%   dvec        the vector of expansion coefficients corresponding to the 
%               bestbasis
%   BS          a BasisSpec object which specifies the best-basis
%   trans       specifies which transform was used for that portion of the
%               signal: 
%                   00 = HGLET with L
%                   01 = HGLET with Lrw
%                   10 = HGLET with Lsym
%                   11 = GHWT
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% specify transform codes
transHsym = [ true,false];
transG    = [ true, true];
transHrw  = [false, true];
transH    = [false,false];

% the cost functional to be used
if ~exist('costfun','var')
    costfun = [];
end
costfun = cost_functional(costfun);

% constants and dmatrix cleanup
if isarray(dmatrixHsym)
    [N,jmax,fcols] = size(dmatrixHsym);
    dmatrixHsym( abs(dmatrixHsym) < 10^2*eps ) = 0;
elseif isarray(dmatrixG)
    [N,jmax,fcols] = size(dmatrixG);
    dmatrixG( abs(dmatrixG) < 10^2*eps ) = 0;
elseif isarray(dmatrixHrw)
    [N,jmax,fcols] = size(dmatrixHrw);
    dmatrixHrw( abs(dmatrixHrw) < 10^2*eps ) = 0;
elseif isarray(dmatrixH)
    [N,jmax,fcols] = size(dmatrixH);
    dmatrixH( abs(dmatrixH) < 10^2*eps ) = 0;
else
    fprintf('\n\nNo coefficient matrices provided.  Exiting now.\n\n');
    return
end
    
% "flatten" dmatrix
if fcols > 1
    if ~exist('flatten','var')
        flatten = 1;
    end
    if isarray(dmatrixHsym)
        dmatrix0Hsym = dmatrixHsym;
        dmatrixHsym = dmatrix_flatten(dmatrixHsym,flatten);
    end
    if isarray(dmatrixG)
        dmatrix0G = dmatrixG;
        dmatrixG = dmatrix_flatten(dmatrixG,flatten);
    end
    if isarray(dmatrixHrw)
        dmatrix0Hrw = dmatrixHrw;
        dmatrixHrw = dmatrix_flatten(dmatrixHrw,flatten);
    end
    if isarray(dmatrixH)
        dmatrix0H = dmatrixH;
        dmatrixH = dmatrix_flatten(dmatrixH,flatten);
    end
end


%% Find the HGLET/GHWT best-basis

% allocate/initialize ==> order matters here
if isarray(dmatrixHsym)
    dvec = dmatrixHsym(:,jmax);
    trans = repmat(transHsym,N,1);
end
if isarray(dmatrixG)
    dvec = dmatrixG(:,jmax);
    trans = repmat(transG,N,1);
end
if isarray(dmatrixHrw)
    dvec = dmatrixHrw(:,jmax);
    trans = repmat(transHrw,N,1);
end
if isarray(dmatrixH)
    dvec = dmatrixH(:,jmax);
    trans = repmat(transH,N,1);
end
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
            
        %%%%% compute the cost of the HGLET-Lsym coefficients
        if isarray(dmatrixHsym)
            costNEW = costfun( dmatrixHsym(indr,j) );
            % change the best basis if the new cost is less expensive
            if costBB >= costNEW - tol
                [costBB, dvec(indr), levlist(indr), trans(indr,:)] = BBchange(costNEW,dmatrixHsym(indr,j),j,transHsym);
            end
        end
        
        %%%%% compute the cost of the GHWT coefficients
        if isarray(dmatrixG)
            costNEW = costfun( dmatrixG(indr,j) );
            % change the best basis if the new cost is less expensive
            if costBB >= costNEW - tol
                [costBB, dvec(indr), levlist(indr), trans(indr,:)] = BBchange(costNEW,dmatrixG(indr,j),j,transG);
            end
        end
        
        %%%%% compute the cost of the HGLET-Lrw coefficients
        if isarray(dmatrixHrw)
            costNEW = costfun( dmatrixHrw(indr,j) );
            % change the best basis if the new cost is less expensive
            if costBB >= costNEW - tol
                [costBB, dvec(indr), levlist(indr), trans(indr,:)] = BBchange(costNEW,dmatrixHrw(indr,j),j,transHrw);
            end
        end
        
        %%%%% compute the cost of the HGLET-L coefficients
        if isarray(dmatrixH)
            costNEW = costfun( dmatrixH(indr,j) );
            % change the best basis if the new cost is less expensive
            if costBB >= costNEW - tol
                [~, dvec(indr), levlist(indr), trans(indr,:)] = BBchange(costNEW,dmatrixH(indr,j),j,transH);
            end
        end
    end
end

transfull = trans;
trans( levlist==0,: ) = [];
levlist( levlist==0 ) = [];

BS = BasisSpec(levlist,[],true,'HGLET-GHWT Best Basis');
BS = levlist2levlengths(GP,BS);


% if we flattened dmatrix, then "unflatten" the expansion coefficients
if fcols > 1
    % create vectors of coefficients (which are zero if the transform's coefficients were not included as function inputs)
    if isarray(dmatrixH)
        dvecH = dmatrix2dvec(dmatrix0H,GP,BS);
    else
        dvecH = zeros(N,fcols);
    end
    if isarray(dmatrixHrw)
        dvecHrw = dmatrix2dvec(dmatrix0Hrw,GP,BS);
    else
        dvecHrw = zeros(N,fcols);
    end
    if isarray(dmatrixHsym)
        dvecHsym = dmatrix2dvec(dmatrix0Hsym,GP,BS);
    else
        dvecHsym = zeros(N,fcols);
    end
    if isarray(dmatrixG)
        dvecG = dmatrix2dvec(dmatrix0G,GP,BS);
    else
        dvecG = zeros(N,fcols);
    end
    
    dvec = bsxfun(@times, dvecHsym,  transfull(:,1) .* ~transfull(:,2)) ...
            + bsxfun(@times, dvecG,  transfull(:,1) .*  transfull(:,2)) ...
          + bsxfun(@times, dvecHrw, ~transfull(:,1) .*  transfull(:,2)) ...
            + bsxfun(@times, dvecH, ~transfull(:,1) .* ~transfull(:,2));
end


end




function [costBB, dvec, levlist, trans] = BBchange(costNEW, dvec, j, trans)
% Change to the new best basis

costBB = costNEW;

n = length(dvec);

levlist = zeros(n,1,'uint8');
levlist(1) = j;

trans = repmat(trans,n,1);
end