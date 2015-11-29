function [dvec,BS,trans,cost,q,T] = HGLET_GHWT_Path_BestBasis_MDL(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GP)
% For a signal on a path graph, select the best basis from several matrices
% of expansion coefficients using the MDL as the cost functional
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
%   q           the quantization precision is delta=2^-q
%   T           the quantization threshold (i.e., quantized expansion
%               coefficients with l~=0 that are <=T will be set to 0)
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

% constants and dmatrix cleanup; prepare to normalize the coefficients
if isarray(dmatrixHsym)
    [N,jmax,dnorm,dmatrixHsym,dmatrixHsym0] = process_array(dmatrixHsym);
else
    dmatrixHsym0 = 0;
end
if isarray(dmatrixG)
    [N,jmax,dnorm,dmatrixG,dmatrixG0] = process_array(dmatrixG);
else
    dmatrixG0 = 0;
end
if isarray(dmatrixHrw)
    [N,jmax,dnorm,dmatrixHrw,dmatrixHrw0] = process_array(dmatrixHrw);
else
    dmatrixHrw0 = 0;
end
if isarray(dmatrixH)
    [N,jmax,dnorm,dmatrixH,dmatrixH0] = process_array(dmatrixH);
else
    dmatrixH0 = 0;
end
if ~exist('N','var')
    fprintf('\n\nNo coefficient matrices provided.  Exiting now.\n\n');
    return
end
    
% extract data
GP = GHWT_Info(GP);
[~,rs,tag] = ExtractData(GP);


%%% MDL costs

% compute the number of bits needed to store 'trans' entries
trans_cost = ceil(log2( isarray(dmatrixHsym) + isarray(dmatrixG) ...
    + isarray(dmatrixHrw) + isarray(dmatrixH)));

% compute the number of bits needed to store levlist entries
levlist_cost = ceil(log2(jmax));

qmin = 1;
qmax = 5;

% define the cost functional
costfun = @(xQ,xQO,xO,q,numregions) MDL_Path(xQ,xQO,xO,q) + numregions*(levlist_cost + trans_cost);


%% Find the HGLET/GHWT best-basis

% set the tolerance for the best-basis search
tol = 10^4*eps;

for q = qmin:qmax
    
    % the precision
    delta = 2^(-q);
    
    % the precision-scaled coefficients and rounded integer coefficients
    [dmatrixHsym,dmatrixHsymQ] = scale_and_round(dmatrixHsym0,dnorm,delta);
    [dmatrixG,dmatrixGQ] = scale_and_round(dmatrixG0,dnorm,delta);
    [dmatrixHrw,dmatrixHrwQ] = scale_and_round(dmatrixHrw0,dnorm,delta);
    [dmatrixH,dmatrixHQ] = scale_and_round(dmatrixH0,dnorm,delta);
    
    for T = 0:2
        % 1) threshold the coefficients
        % 2) initialize the precision-scaled best basis coefficcients
        % 3) initialize the quantized (integer) best basis coefficients
        % 4) initialize the orthonormal quantized best basis coefficients
        % 5) initialize the transform specifications        
        if isarray(dmatrixHsym)
            dmatrixHsymQ = dmatrix_threshold(dmatrixHsymQ,T,tag);
            [dvecO,dvecQ,dvecQO,trans] = initialize_bestbasis(dmatrixHsym,dmatrixHsymQ,transHsym);
        end
        if isarray(dmatrixG)
            dmatrixGQ = dmatrix_threshold(dmatrixGQ,T,tag);
            [dvecO,dvecQ,dvecQO,trans] = initialize_bestbasis(dmatrixG,dmatrixGQ,transG);
        end
        if isarray(dmatrixHrw)
            dmatrixHrwQ = dmatrix_threshold(dmatrixHrwQ,T,tag);
            [dvecO,dvecQ,dvecQO,trans] = initialize_bestbasis(dmatrixHrw,dmatrixHrwQ,transHrw);
        end
        if isarray(dmatrixH)
            dmatrixHQ = dmatrix_threshold(dmatrixHQ,T,tag);
            [dvecO,dvecQ,dvecQO,trans] = initialize_bestbasis(dmatrixH,dmatrixHQ,transH);
        end
        
        % the levels_list description of the best basis
        levlist = jmax*ones(N,1,'uint8');
        
        % perform the basis search
        for j = jmax-1:-1:1
            regioncount = nnz(rs(:,j))-1;
            for r = 1:regioncount
                indr = rs(r,j):rs(r+1,j)-1;
                if ~isempty(indr)
                    % compute the cost of the current best basis
                    costBB = costfun(dvecQ(indr),dvecQO(indr),dvecO(indr),q,nnz(levlist(indr)));

                    %%%%% compute the cost of the HGLET-Lsym coefficients
                    if isarray(dmatrixHsym)
                        costNEW = costfun(dmatrixHsymQ(indr,j),dmatrixHsymQ(indr,j),dmatrixHsym(indr,j),q,1);

                        % change the best basis if the new cost is less expensive
                        if costBB >= costNEW - tol
                            [costBB, dvecQ(indr), dvecQO(indr), dvecO(indr), levlist(indr), trans(indr,:)] = BBchange(costNEW,dmatrixHsymQ(indr,j),dmatrixHsymQ(indr,j),dmatrixHsym(indr,j),j,transHsym);
                        end
                    end

                    %%%%% compute the cost of the GHWT coefficients
                    if isarray(dmatrixG)
                        costNEW = costfun(dmatrixGQ(indr,j),dmatrixGQ(indr,j),dmatrixG(indr,j),q,1);

                        % change the best basis if the new cost is less expensive
                        if costBB >= costNEW - tol
                            [costBB, dvecQ(indr), dvecQO(indr), dvecO(indr), levlist(indr), trans(indr,:)] = BBchange(costNEW,dmatrixGQ(indr,j),dmatrixGQ(indr,j),dmatrixG(indr,j),j,transG);
                        end
                    end

                    %%%%% compute the cost of the HGLET-Lrw coefficients
                    if isarray(dmatrixHrw)
                        [dQO,dO] = Lrw_orth_coeffs(dmatrixHrwQ(indr,j),dmatrixHrw(indr,j));
                        costNEW = costfun(dmatrixHrwQ(indr,j),dQO,dO,q,1);

                        % change the best basis if the new cost is less expensive
                        if costBB >= costNEW - tol
                            [costBB, dvecQ(indr), dvecQO(indr), dvecO(indr), levlist(indr), trans(indr,:)] = BBchange(costNEW,dmatrixHrwQ(indr,j),dQO,dO,j,transHrw);
                        end
                    end

                    %%%%% compute the cost of the HGLET-L coefficients
                    if isarray(dmatrixH)
                        costNEW = costfun(dmatrixHQ(indr,j),dmatrixHQ(indr,j),dmatrixH(indr,j),q,1);

                        % change the best basis if the new cost is less expensive
                        if costBB >= costNEW - tol
                            [~, dvecQ(indr), dvecQO(indr), dvecO(indr), levlist(indr), trans(indr,:)] = BBchange(costNEW,dmatrixHQ(indr,j),dmatrixHQ(indr,j),dmatrixH(indr,j),j,transH);
                        end
                    end
                end
            end
        end

        % compare the cost for this q and T versus the current best basis
        cost = costfun(dvecQ,dvecQO,dvecO,q,nnz(levlist))+q;
        if ~exist('cost0','var') || cost < cost0
            cost0  = cost;
            q0 = q;
            T0 = T;
            trans0 = trans;
            levlist0 = levlist;
        end
    end % T for-loop
end % q for-loop

% specify the best basis
cost = cost0;
trans = trans0;
levlist = levlist0;
trans( levlist==0,: ) = [];
levlist( levlist==0 ) = [];
BS = BasisSpec(levlist,[],true,'HGLET-GHWT Best Basis');
BS = levlist2levlengths(GP,BS);

% return quantized & thresholded coefficients ==> use the best q and T
delta = 2^(-q0);
dmatrixHsym = scale_round_threshold_rescale(dmatrixHsym0,dnorm,delta,T0,tag);
dmatrixG = scale_round_threshold_rescale(dmatrixG0,dnorm,delta,T0,tag);
dmatrixHrw = scale_round_threshold_rescale(dmatrixHrw0,dnorm,delta,T0,tag);
dmatrixH = scale_round_threshold_rescale(dmatrixH0,dnorm,delta,T0,tag);

dvec = dmatrices2dvec(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GP,BS,trans);


end




function [cost, dvecQ, dvecQO, dvecO, levlist, trans] = BBchange(cost,dvecQ,dvecQO,dvecO,j,trans)
% Change to the new best basis

n = length(dvecQ);
levlist = zeros(n,1,'uint8');
levlist(1) = j;
trans = repmat(trans,n,1);
end




function [dQO,dO] = Lrw_orth_coeffs(dQ,d)
% Synthesize the signal using Lrw coefficients and quantized Lrw
% coefficients

n = length(dQ);
D = 2*ones(n,1);
D([1,n]) = 1;

dQO = bsxfun(@times, D.^(-0.5), idct1(dQ));
dO  = bsxfun(@times, D.^(-0.5), idct1(d));
end




function [N,jmax,dnorm,dmatrix,dmatrix0] = process_array(dmatrix)
% constants and dmatrix cleanup; prepare to normalize the coefficients

[N,jmax,~] = size(dmatrix);
dmatrix( abs(dmatrix) < 10^2*eps ) = 0;
dnorm = norm(dmatrix(:,end),2)/sqrt(N);
dmatrix0 = dmatrix;
end




function [dmatrix,dmatrixQ] = scale_and_round(dmatrix0,dnorm,delta)
% scale the coefficients to precision delta and round them

dmatrix  = dmatrix0/dnorm/2/delta;
dmatrixQ = round(dmatrix); 
end




function dmatrixQ = dmatrix_threshold(dmatrixQ,T,tag)
% threshold the quantized coefficients

dmatrixQ( -T <= dmatrixQ & dmatrixQ <= T & tag ~= 0 ) = 0;
end




function [dvecO,dvecQ,dvecQO,trans] = initialize_bestbasis(dmatrix,dmatrixQ,trans)
% initialize the best basis to the finest level of the current dictionary

[N,~] = size(dmatrix);
dvecO = dmatrix(:,end);
dvecQ = dmatrixQ(:,end);
dvecQO = dmatrixQ(:,end);
trans = repmat(trans,N,1);
end




function dmatrix = scale_round_threshold_rescale(dmatrix0,dnorm,delta,T,tag)
% scale, round, threshold, and then rescale the coefficients

dmatrix = round(dmatrix0/dnorm/2/delta);
if isarray(dmatrix)
    dmatrix( -T <= dmatrix & dmatrix <= T & tag ~= 0 ) = 0;
    dmatrix = dmatrix*dnorm*2*delta;
end
end