function fig = GraphSig_Plot_BlackWhite(G,size0,size1)
% Display a GraphSig object so that it looks good in color and in black &
% white
%
% Input
%   G       a GraphSig object
%   size0   the smallest node size
%   size1   the largest node size
%
% Output
%    fig    a plot figure
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% make sure the graph is 2-D or 3-D
if G.dim ~= 2 && G.dim ~= 3
    fig = [];
    fprintf('\n\nThe graph must be 2-D or 3-D.  Exiting now.\n\n');
    return
end

% plot the graph
fig = figure;
gplot3(G.W,G.xy,'LineStyle','-','Color','k');
hold on

% find the biggest value (in magnitude)
fmax = max(abs(G.f));

% the vector of node sizes
if ~exist('size0','var')
    size0 = 10;
end
if ~exist('size1','var')
    size1 = 50;
end
sizes = size0 + size1*abs(G.f)/fmax;

% only plot the nonzero nodes
nz = find(G.f);

% nodes to be encircled
ec = find(G.f > 0);

% plot the nodes
if G.dim == 2
    scatter(G.xy(nz,1),G.xy(nz,2),sizes(nz),G.f(nz),'o','filled','MarkerEdgeColor','none');
    scatter(G.xy(ec,1),G.xy(ec,2),sizes(ec),G.f(ec),'o','filled','MarkerEdgeColor','k');
else
    scatter3(G.xy(nz,1),G.xy(nz,2),G.xy(nz,3),sizes(nz),G.f(nz),'o','filled','MarkerEdgeColor','none');
    scatter3(G.xy(ec,1),G.xy(ec,2),G.xy(ec,3),sizes(ec),G.f(ec),'o','filled','MarkerEdgeColor','k');
end

% the colormap
cmap = ones(64,3);
cmap(1:22,1:3)  = [(0:21)'/22, zeros(22,2)];    
cmap(22:43,2:3) = [(0:21)'/22, zeros(22,1)];
cmap(43:64,3) = (0:21)'/22;
colormap(cmap);
set(gca, 'CLim', [-fmax, fmax]);


set(gcf,'color','w');

set(gca, 'FontSize',10);



end