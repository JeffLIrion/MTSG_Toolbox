% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



load('MN_MutGauss.mat','G');
G = EditPlotSpecs(G,'symm verbatim{{axis equal;}}');

% convert to inverse Euclidean edge weights
G = Adj2InvEuc(G);

% partition the graph with the Fiedler vectors of Lrw
GP = PartitionTreeFiedler(G,1);

% extract RegionStarts data
[~,rs] = ExtractData(GP);

% Figure 3a
G3a = HGLET_BasisVector(G,GP,6,2);
[j1,k1,l1] = GHWT_jkl(GP,6,2);

% Figure 3b
G3b = HGLET_BasisVector(G,GP,rs(4,4)+1,4);
[j2,k2,l2] = GHWT_jkl(GP,rs(4,4)+1,4);


%%% display the figures
figure;

% the colormap to be used
cmap = ones(64,3);
cmap(1:22,1:3)  = [(0:21)'/22, zeros(22,2)];    
cmap(22:43,2:3) = [(0:21)'/22, zeros(22,1)];
cmap(43:64,3) = (0:21)'/22;
colormap(cmap);

% Figure 3a
s1 = subplot(1,2,1);
axis equal;
axis off;
cmax = max(abs(G3a));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\phi^{%d}_{%d,%d}$$',j1,k1,l1),'interpreter','latex','FontSize',16);
GraphSig_Plot_BlackWhite(G3a);

% grab the plot and put it in the figure
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% Figure 3b
s1 = subplot(1,2,2);
axis equal;
axis off;
cmax = max(abs(G3b));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\phi^{%d}_{%d,%d}$$',j2,k2,l2),'interpreter','latex','FontSize',16);
GraphSig_Plot_BlackWhite(G3b);

% grab the plot and put it in the figure
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


set(gcf,'color','w');