function [pm,v] = PartitionFiedlerSweep(W,~,v)
% Partition the vertices of a graph by sweeping through the entries of the 
% Fiedler vector.  
%
% Inputs
%   W           the edge weight matrix
%   ~           if a 2nd input argument is given, use L_rw
%   v           the Fiedler vector
%
% Outputs
%   pm          a vector of 1's and -1's
%   v           the Fiedler vector
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% Easy case: the Fiedler vector is provided

if exist('v','var')
    % partition using the Fiedler vector and troubleshoot any potential issues
    [pm,v] = PartitionFiedlerSweep_pm(v);
    
    if length(W) == length(v)
        val = v'*(diag(sum(W))-W)*v;
        pm = PartitionFiedlerSweep_troubleshooting(pm,v,W,val);
    else
        pm = PartitionFiedlerSweep_troubleshooting(pm);
    end
    
    return
end


%% Preliminaries

N = length(W);
eigs_flag = 0;

sigma = eps;
cutoff = 128; % this value could be changed...


%% handle the case when there are 2 nodes
if N == 2
    pm = [1;-1];
    v = pm/sqrt(2);
    return
end


%% L (unnormalized)
if nargin == 1 || min(sum(W)) < 10^3*eps
    if N > cutoff
        opts.issym = 1;
        opts.v0 = ones(N,1)/sqrt(N);
        warning('error','MATLAB:nearlySingularMatrix');
        warning('error','MATLAB:eigs:SigmaNearExactEig');
        
        try
            [v,val,eigs_flag] = eigs(diag(sum(W))-W,2,sigma,opts);
        catch
            eigs_flag = 2;
        end
        
        if eigs_flag == 0
            [~,ind] = max(diag(val));
            v = v(:,ind);
            val = val(ind,ind);
        end
    end
    if N <= cutoff || eigs_flag ~= 0
        %%% Fiedler vector via svd
        [v,val,~] = svd(full(diag(sum(W))-W));
        v = v(:,end-1);
        val = val(end-1,end-1);
    end
    
    
%% L_rw
elseif nargin == 2
    if N > cutoff
        opts.issym = 1;
        opts.v0 = ones(N,1)/sqrt(N);
        warning('error','MATLAB:nearlySingularMatrix');
        warning('error','MATLAB:eigs:SigmaNearExactEig');
        
        try
            [v,val,eigs_flag] = eigs(diag(sum(W))-W,diag(sum(W)),2,sigma,opts);
        catch
            eigs_flag = 2;
        end
        
        if eigs_flag == 0
            [~,ind] = max(diag(val));
            v = v(:,ind);
            v = (full(sum(W,2)).^(-0.5)) .* v;
            val = val(ind,ind);
        end
    end
    if N <= cutoff || eigs_flag ~= 0
        %%% Fiedler vector via svd
        [v,val,~] = svd( full( bsxfun(@times, bsxfun(@times, full(sum(W,2)).^(-0.5), diag(sum(W)) - W), (full(sum(W,1))).^(-0.5) ) ) );
        v = v(:,end-1);
        v = (full(sum(W,2)).^(-0.5)) .* v;
        val = val(end-1,end-1);
    end
end


% partition using the Fiedler vector and troubleshoot any potential issues
[pm,v] = PartitionFiedlerSweep_pm(v);
pm = PartitionFiedlerSweep_troubleshooting(pm,v,W,val);


end




function [pm,v] = PartitionFiedlerSweep_pm(v)
% Partition based on the Fiedler vector 'v'

% indices corresponding to 25% and 75% of the vector length
N = length(v);
ind25 = ceil(0.25*N);
ind75 = floor(0.75*N);

% find the maximum difference between sorted entries of the Fiedler vector
x = sort(v,'ascend');
y = diff(x);
[~,IX] = max(y(ind25:ind75));
IX = IX+ind25-1;

% the cutoff
cutoff = (x(IX)+x(IX+1))/2;
if cutoff <= min(v) || cutoff >= max(v)
    cutoff = mean(v);
end

% create an assignment vector of +1s and -1s
pm = 2*(v >= cutoff)-1;

% make sure the first point is assigned to region 1 (not -1)
if pm(1) < 0
    pm = -pm;
end
end



function pm = PartitionFiedlerSweep_troubleshooting(pm,v,W,val)
% Troubleshoot potential issues with the partitioning

N = length(pm);
tol = 10^3*eps;

if sum(pm < 0) == 0 || sum(pm > 0) == 0 || sum(abs(pm)) < N
    if nargin == 4
        % Case 1: it is not connected
        if val < tol
            pm = 2*(abs(v) > tol) - 1;
            while sum(pm < 0) == 0 || sum(pm > 0) == 0
                tol = 10*tol;
                pm = 2*(abs(v) > tol) - 1;
                if tol > 1
                    pm(1:ceil(N/2)) = 1;
                    pm(ceil(N/2)+1:N) = -1;
                end
            end

        % Case 2: it is connected
        else
            pm = 2*(v >= mean(v)) - 1;

            if sum(abs(pm)) < N
                % assign the near-zero points based on the values of v at their
                % neighbor nodes
                pm0 = (1:N)';
                pm0 = pm0(pm == 0);
                pm(pm0) = (W(pm0,:)*v > tol) - (W(pm0,:)*v < -tol);

                % assign any remaining zeros to the group with fewer members
                pm(pm == 0) = (sum(pm > 0) - sum(pm < 0)) >= 0;
            end
        end
    end

    % if one region has no points
    if sum(pm < 0) == 0 || sum(pm > 0) == 0
        pm(1:ceil(N/2)) = 1;
        pm(ceil(N/2)+1:N) = -1;
    end    
end

% make sure that the first point is assigned as a 1
if pm(1) < 0
    pm = -pm;
end
end