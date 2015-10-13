function fig = Display_GHWT_Basis_Vectors_c2f(G,GP)
% Display a figure of the GHWT basis vectors in their coarse-to-fine
% arrangement
%
% Input
%   G       a GraphSig object
%   GP      a GraphPart object
%
% Output
%   fig     a figure of the GHWT basis vectors (coarse-to-fine dictionary)
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


% generate the GHWT coarse-to-fine dictionary
figure;
for j = 1:jmax
    for row = 1:N
        
        % generate the basis vector
        [~,k,l] = GHWT_jkl(GP,row,j);
        dmatrix = zeros(N,jmax);
        dmatrix(row,j) = 1;
        [dvec,BS] = dmatrix2dvec(dmatrix,GP);
        [~,Gout] = GHWT_Synthesis(dvec,GP,BS,G);
        
        % 1-D case: specify the color of the nodes
        if dim == 1
            if l == 0
                Gout = EditPlotSpecs(Gout,'linecolork',1);
            elseif l == 1
                Gout = EditPlotSpecs(Gout,'linecolorr',1);
            else
                Gout = EditPlotSpecs(Gout,'linecolorb',1);
            end
        end
        
        % plot the basis vector
        s1 = subplot(jmax,N,N*(j-1)+row);
        axis equal;
        axis off;
        if dim == 1
            axis([0.9, N+0.1, -1, 1]);
        else
            set(gca, 'CLim', [-1, 1]);
        end
        title(sprintf('$$\\psi^{%d}_{%d,%d}$$',j-1,k,l),'interpreter','latex','FontSize',16);
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