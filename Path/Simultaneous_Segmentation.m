function [G,GP,BS,trans,q,T] = Simultaneous_Segmentation(G,Nmin,GP,BS)
% Use the MDL as the best basis cost functional in order to simultaneously
% segment and denoise a 1-D signal on a path graph.  
%
% Input
%   G           a GraphSig object
%   Nmin        the minimum segment length in the denoised signal
%               (default = [N/50])
%   GP          initial recursive graph partitioning (optional)
%   BS          initial basis (optional)
%
% Output
%   G           the denoised signal
%   GP          the recursive graph partitioning
%   BS          the basis for the denoised signal
%   trans       specifies which transform was used for that portion of the
%               signal: 
%                   00 = HGLET with L    ==> blue
%                   01 = HGLET with Lrw  ==> red
%                   10 = HGLET with Lsym ==> black
%   q           the quantization precision is delta=2^-q
%   T           the quantization threshold (i.e., quantized expansion
%               coefficients with l~=0 that are <=T will be set to 0)
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% initialize variables
iter = 0;
cost0 = 0;
samecost = 0;

% partition the graph
if exist('GP','var') && exist('BS','var')
    GP = PartialTreePathNCut(G,GP,BS,'choose');
else
    GP = PartitionTreePathNCut(G);
end

while samecost < 3
    iter = iter+1;

    % re-partition
    if iter ~= 1
        GP1 = GP;
        if samecost == 0
            GP = PartialTreePathNCut(G,GP1,BS,'choose');
        elseif samecost == 1
            GP = PartialTreePathNCut(G,GP1,BS,'yes');
        elseif samecost == 2
            GP = PartialTreePathNCut(G,GP1,BS,'no');
        end
    end

    % analyze the signal
    [dH,dHrw,dHsym] = HGLET_Path_Analysis_All(G,GP);

    % find the best basis
    [dvec,BS,trans,cost,q,T] = HGLET_GHWT_Path_BestBasis_MDL(dH,dHrw,dHsym,0,GP);

    % cost INCREASED ==> revert to previous state
    if iter > 1 && cost > cost0 + 10^-3
        fprintf('\n\nCOST INCREASED!!!!! ==> revert to previous state\n\n');
        G = G0;
        GP = GP0;
        BS = BS0;
        trans = trans0;
        q = q0;
        T = T0;
        samecost = samecost+1;

    else
        if iter > 1
            % cost stayed the same
            if abs(cost-cost0) <= 10^-3
                samecost = samecost+1;

            % cost decreased
            else
                samecost = 0;
            end
        end

        % cost didn't increase 
        %%% 1) bacqup the current state
        G0 = G;
        GP0 = GP;
        BS0 = BS;
        trans0 = trans;
        cost0 = cost;
        q0 = q;
        T0 = T;

        %%% 2) modify the edge weights
        G = Modify_Edge_Weights(G,GP,BS);
    end

    % visualize the best basis
    if iter > 1
        close(gcf);
    end
    Path_BasisVisual(G,GP,BS,trans);
    title(sprintf('Round %d (samecost=%d)',iter,samecost));
    
    % plot the modified edge weights on top of the best basis visualization
    W = ExtractData(G);
    W = full(diag(W,1)-1)*norm(G,Inf)/4+mean(ylim);
    x = 0.5+(1:length(W))';
    plot(x,W,'-g');
    pause(0.1);

    % print the cost for the current round
    fprintf('\n\nRound %2d cost: %10f  (q=%d, T=%d)\n',iter,cost,q,T);
end

fprintf('\n\nFinished after %d iterations.\n\n',iter)

% synthesize from the MDL-guided best basis coefficients
[~,G] = HGLET_GHWT_Path_Synthesis(dvec,GP,BS,trans,G);

% plot the denoised and segmented signal
close(gcf);
if exist('Nmin','var')
    Path_BasisVisual(G,GP,BS,trans,Nmin);
else
    Path_BasisVisual(G,GP,BS,trans,round(length(G)/50));
end


end




function G = Modify_Edge_Weights(G,GP,BS)
% For every edge that is cut in the current best basis, cut two edges that
% are 5% and 10% to the left and right of it.  (The percentages are in
% reference to the previous/following cut edge.)

% extract data
[~,xy,f,name,plotspecs] = ExtractData(G);
N = length(G);
[~,levlengths] = ExtractData(BS,GP);

% unweighted path graph adjacency matrix
W = ExtractData(Gpath(N));

% generate the regionstarts for the basis at hand
rs = [0; cumsum(double(levlengths))]+1;

for row = 2:length(rs)-1
    % the starting nodes of the regions
    rs1 = rs(row-1);
    rs2 = rs(row);
    rs3 = rs(row+1);

    % the 10% edges to be cut
    n1 = max([1, round((rs2-rs1-1)/10)]);
    n2 = max([1, round((rs3-rs2-1)/10)]);

    % W(rs2-1,rs2) is the edge that is cut in the current best basis
    W(rs2-n1-1,rs2-n1) = 0;
    W(rs2-n1,rs2-n1-1) = 0;
    W(rs2+n2-1,rs2+n2) = 0;
    W(rs2+n2,rs2+n2-1) = 0;

    % the 5% edges to be cut
    n1 = max([1, round((rs2-rs1-1)/20)]);
    n2 = max([1, round((rs3-rs2-1)/20)]);

    % W(rs2-1,rs2) is the edge that is cut in the current best basis
    W(rs2-n1-1,rs2-n1) = 0;
    W(rs2-n1,rs2-n1-1) = 0;
    W(rs2+n2-1,rs2+n2) = 0;
    W(rs2+n2,rs2+n2-1) = 0;
end

% rebuild G
G = GraphSig(W,xy,f,name,plotspecs);
end