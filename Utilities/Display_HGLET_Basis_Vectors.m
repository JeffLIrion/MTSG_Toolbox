function fig = Display_HGLET_Basis_Vectors(G,GP)
% Display a figure of the HGLET basis vectors
%
% Input
%   G       a GraphSig object
%   GP      a GraphPart object
%
% Output
%   fig     a figure of the HGLET basis vectors
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% partition the graph, if necessary
if ~exist('GP','var')
    GP = PartitionTreeFiedler(G);
end

% constants
[~,xy] = ExtractData(G);
[~,dim] = size(xy);
[~,rs] = ExtractData(GP);
[N,jmax] = size(rs);
N = N-1;

if dim == 1
    G = EditPlotSpecs(G,'size10 linewidth2 stem');
end

% make sure the graph is small enough
if N > 10
    fprintf('\n\nThis graph is too big to display the full dictionary.  Exiting now.\n\n');
    fig = [];
    return
end


% generate the HGLET dictionary
fig = figure;
for j = 1:jmax
    for row = 1:N
        
        % generate the basis vector
        [~,k,l] = HGLET_jkl(GP,row,j);
        dmatrix = zeros(N,jmax);
        dmatrix(row,j) = 1;
        [dvec,BS] = dmatrix2dvec(dmatrix,GP);
        [~,Gout] = HGLET_Synthesis(dvec,GP,BS,G);
        
        % plot the basis vector
        s1 = subplot(jmax,N,N*(j-1)+row);
        axis equal;
        axis off;
        if dim == 1
            axis([0.9, N+0.1, -1, 1]);
        else
            set(gca, 'CLim', [-1, 1]);
        end
        title(sprintf('$$\\phi^{%d}_{%d,%d}$$',j-1,k,l),'interpreter','latex','FontSize',16);
        GraphSig_Plot(Gout);
        
        % grab the plot and put it in the figure with all the basis vectors
        ax1 = gca;
        h1 = gcf;
        fig1 = get(ax1,'children');
        copyobj(fig1,s1);
        close(h1);
    end
end

set(gcf,'color','w');


end