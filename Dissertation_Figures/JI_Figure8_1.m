% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



figure;
colormap('jet');

load('MN_MutGauss.mat');
G0 = G;
load('MN_MutGauss_SNR5.mat');

% Figure 8.1a
s1 = subplot(1,2,1);
axis equal;
xlim([-97.5, -89]);
ylim([42, 50]);
set(gca, 'CLim', [min(G), max(G)]);
GraphSig_Plot(G0);

% grab the plot and put it in the figure
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% Figure 8.1b
s1 = subplot(1,2,2);
axis equal;
xlim([-97.5, -89]);
ylim([42, 50]);
set(gca, 'CLim', [min(G), max(G)]);
GraphSig_Plot(G);

% grab the plot and put it in the figure
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


set(gcf,'color','w');