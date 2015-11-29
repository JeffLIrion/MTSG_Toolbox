% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% load the graph
load('MN_MutGauss','G');
G = EditPlotSpecs(G,'symm notitle nocolorbar sortnodes verbatim{{axis equal; axis off;}}');

% partition the graph
GP = PartitionTreeFiedler(G,1);
N = length(G);
[~,rs] = ExtractData(GP);


%% level j=0 (global) basis vectors
BS = LevelBasisSpec(GP,0);
drow = [2,11,21,51];
dvec = zeros(N,4);
dvec(drow(1),1) = 1;
dvec(drow(2),2) = 1;
dvec(drow(3),3) = 1;
dvec(drow(4),4) = 1;
fH = HGLET_Synthesis(dvec,GP,BS,G);

[j,k,l] = HGLET_jkl(GP,drow(1),1); 
H1 = EditName(ReplaceData(G,fH(:,1)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = HGLET_jkl(GP,drow(2),1); 
H2 = EditName(ReplaceData(G,fH(:,2)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = HGLET_jkl(GP,drow(3),1); 
H3 = EditName(ReplaceData(G,fH(:,3)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = HGLET_jkl(GP,drow(4),1); 
H4 = EditName(ReplaceData(G,fH(:,4)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));


%% level j=1 basis vectors
BS = LevelBasisSpec(GP,1);
drow = [2,11,rs(2,2)+1,rs(2,2)+15];
dvec = zeros(N,4);
dvec(drow(1),1) = 1;
dvec(drow(2),2) = 1;
dvec(drow(3),3) = 1;
dvec(drow(4),4) = 1;
fH = HGLET_Synthesis(dvec,GP,BS,G);

[j,k,l] = HGLET_jkl(GP,drow(1),2); 
H5 = EditName(ReplaceData(G,fH(:,1)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = HGLET_jkl(GP,drow(2),2); 
H6 = EditName(ReplaceData(G,fH(:,2)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = HGLET_jkl(GP,drow(3),2); 
H7 = EditName(ReplaceData(G,fH(:,3)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = HGLET_jkl(GP,drow(4),2); 
H8 = EditName(ReplaceData(G,fH(:,4)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));


%% level j=2 basis vectors
BS = LevelBasisSpec(GP,2);
drow = [2,rs(2,3)+2,rs(3,3)+3,rs(4,3)+4];
dvec = zeros(N,4);
dvec(drow(1),1) = 1;
dvec(drow(2),2) = 1;
dvec(drow(3),3) = 1;
dvec(drow(4),4) = 1;
fH = HGLET_Synthesis(dvec,GP,BS,G);

[j,k,l] = HGLET_jkl(GP,drow(1),3); 
H9 = EditName(ReplaceData(G,fH(:,1)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = HGLET_jkl(GP,drow(2),3); 
H10 = EditName(ReplaceData(G,fH(:,2)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = HGLET_jkl(GP,drow(3),3); 
H11 = EditName(ReplaceData(G,fH(:,3)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = HGLET_jkl(GP,drow(4),3); 
H12 = EditName(ReplaceData(G,fH(:,4)),sprintf('$$\\phi^{%d}_{%d,%d}$$',j,k,l));


%% generate HGLET figures in Matlab

h = figure;

s1 = subplot(3,4,1);
axis equal;
axis off;
cmax = max(abs(H1));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H1);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H1);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,2);
axis equal;
axis off;
cmax = max(abs(H2));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H2);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H2);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,3);
axis equal;
axis off;
cmax = max(abs(H3));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H3);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H3);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,4);
axis equal;
axis off;
cmax = max(abs(H4));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H4);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H4);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,5);
axis equal;
axis off;
cmax = max(abs(H5));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H5);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H5);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,6);
axis equal;
axis off;
cmax = max(abs(H6));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H6);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H6);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,7);
axis equal;
axis off;
cmax = max(abs(H7));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H7);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H7);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,8);
axis equal;
axis off;
cmax = max(abs(H8));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H8);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H8);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,9);
axis equal;
axis off;
cmax = max(abs(H9));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H9);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H9);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,10);
axis equal;
axis off;
cmax = max(abs(H10));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H10);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H10);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,11);
axis equal;
axis off;
cmax = max(abs(H11));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H11);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H11);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,12);
axis equal;
axis off;
cmax = max(abs(H12));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(H12);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(H12);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

colormap('jet');
set(gcf,'color','w');


%% display the coarse-to-fine and fine-to-coarse trees

PartitionTreeDisplay(GP);
title('Partitioning Tree');

[c2f,f2c] = HGLET_GHWT_BasisCount(GP);
c2f = sprintf('10^%d = 2^%d',floor(-c2f),floor(-c2f/log10(2)));

fprintf('\n\nThere are (rounded due to machine arithmetic):\n\n');
fprintf('* %s choosable bases from the HGLET (Lrw) dictionary\n\n',c2f);