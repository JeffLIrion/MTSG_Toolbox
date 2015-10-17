% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% NOTE: Due to updates in the code, the basis functions generated herein
% may differ from those in the article.  In particular, due to changes in
% the recursive partitioning code, 
% 1) the signs are flipped in the Northeast corner of Fig. 1b
% 2) Fig. 1c now corresponds to psi(j=2,k=0,l=3) instead of 
%    psi(j=2,k=0,l=2)

load('MN_MutGauss.mat');
G = EditPlotSpecs(G,'symm verbatim{{axis equal;}}');

% convert to inverse Euclidean weights
G = Adj2InvEuc(G);

% partition the graph
GP = PartitionTreeFiedler(G,1);

% extract the regionstarts data
[~,rs] = ExtractData(GP);

% Fig 1a
G1a = GHWT_BasisVector(G,GP,2,1);
[j1,k1,l1] = GHWT_jkl(GP,2,1);

% Fig 1b
G1b = GHWT_BasisVector(G,GP,8,1);
[j2,k2,l2] = GHWT_jkl(GP,8,1);

% Fig 1c
G1c = GHWT_BasisVector(G,GP,4,3);
[j3,k3,l3] = GHWT_jkl(GP,4,3);

% Fig 1d
G1d = GHWT_BasisVector(G,GP,rs(8,4)+3,4);
[j4,k4,l4] = GHWT_jkl(GP,rs(8,4)+3,4);


%%% display the figures

figure

% Figure 1a
s1 = subplot(2,2,1);
axis equal;
axis off;
cmax = max(abs(G1a));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\psi^{%d}_{%d,%d}$$',j1,k1,l1),'interpreter','latex','FontSize',16);
GraphSig_Plot(G1a);

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
cmax = max(abs(G1b));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\psi^{%d}_{%d,%d}$$',j2,k2,l2),'interpreter','latex','FontSize',16);
GraphSig_Plot(G1b);

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
cmax = max(abs(G1c));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\psi^{%d}_{%d,%d}$$',j3,k3,l3),'interpreter','latex','FontSize',16);
GraphSig_Plot(G1c);

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
cmax = max(abs(G1d));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\psi^{%d}_{%d,%d}$$',j4,k4,l4),'interpreter','latex','FontSize',16);
GraphSig_Plot(G1d);

% grab the plot and put it in the figure with all the basis vectors
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


set(gcf,'color','w');