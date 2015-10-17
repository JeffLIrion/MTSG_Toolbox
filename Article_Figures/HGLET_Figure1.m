% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% NOTE: Due to updates in the code, the signs of the basis functions 
% generated herein may differ from those in the article.  

load('MN_MutGauss.mat');
% G = EditPlotSpecs(G,'symm verbatim{{axis equal;}}');

% convert to inverse Euclidean weights
G = Adj2InvEuc(G);

% partition the graph
GP = PartitionTreeFiedler(G,1);

% extract the regionstarts data
[~,rs] = ExtractData(GP);

% Fig 1a
H1a = HGLET_BasisVector(G,GP,rs(2,2)+4,2,1);
[j1,k1,l1] = HGLET_jkl(GP,rs(2,2)+4,2);

% Fig 1b
H1b = HGLET_BasisVector(G,GP,rs(3,3)+1,3,1);
[j2,k2,l2] = HGLET_jkl(GP,rs(3,3)+1,3);

% Fig 1c
H1c = HGLET_BasisVector(G,GP,rs(6,4)+2,4,1);
[j3,k3,l3] = HGLET_jkl(GP,rs(6,4)+2,4);

% Fig 1d
H1d = HGLET_BasisVector(G,GP,rs(6,4)+9,4,1);
[j4,k4,l4] = HGLET_jkl(GP,rs(6,4)+9,4);


%%% display the figures

figure

% Figure 1a
s1 = subplot(2,2,1);
axis equal;
axis off;
cmax = max(abs(H1a));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\phi^{%d}_{%d,%d}$$',j1,k1,l1),'interpreter','latex','FontSize',16);
GraphSig_Plot(H1a);

% grab the plot and put it in the figure with all the basis vectors
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% Figure 1b
s1 = subplot(2,2,2);
axis equal;
axis off;
cmax = max(abs(H1b));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\phi^{%d}_{%d,%d}$$',j2,k2,l2),'interpreter','latex','FontSize',16);
GraphSig_Plot(H1b);

% grab the plot and put it in the figure with all the basis vectors
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% Figure 1c
s1 = subplot(2,2,3);
axis equal;
axis off;
cmax = max(abs(H1c));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\phi^{%d}_{%d,%d}$$',j3,k3,l3),'interpreter','latex','FontSize',16);
GraphSig_Plot(H1c);

% grab the plot and put it in the figure with all the basis vectors
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

% Figure 1d
s1 = subplot(2,2,4);
axis equal;
axis off;
cmax = max(abs(H1d));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\phi^{%d}_{%d,%d}$$',j4,k4,l4),'interpreter','latex','FontSize',16);
GraphSig_Plot(H1d);

% grab the plot and put it in the figure with all the basis vectors
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


set(gcf,'color','w');