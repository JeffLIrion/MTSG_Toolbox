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
fG = GHWT_Synthesis(dvec,GP,BS);

[j,k,l] = GHWT_jkl(GP,drow(1),1);
G1 = EditName(ReplaceData(G,fG(:,1)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = GHWT_jkl(GP,drow(2),1);
G2 = EditName(ReplaceData(G,fG(:,2)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = GHWT_jkl(GP,drow(3),1);
G3 = EditName(ReplaceData(G,fG(:,3)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = GHWT_jkl(GP,drow(4),1);
G4 = EditName(ReplaceData(G,fG(:,4)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));


%% level j=1 basis vectors
BS = LevelBasisSpec(GP,1);
drow = [2,11,rs(2,2)+1,rs(2,2)+15];
dvec = zeros(N,4);
dvec(drow(1),1) = 1;
dvec(drow(2),2) = 1;
dvec(drow(3),3) = 1;
dvec(drow(4),4) = 1;
fG =  GHWT_Synthesis(dvec,GP,BS);

[j,k,l] = GHWT_jkl(GP,drow(1),2);
G5 = EditName(ReplaceData(G,fG(:,1)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = GHWT_jkl(GP,drow(2),2);
G6 = EditName(ReplaceData(G,fG(:,2)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = GHWT_jkl(GP,drow(3),2);
G7 = EditName(ReplaceData(G,fG(:,3)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = GHWT_jkl(GP,drow(4),2);
G8 = EditName(ReplaceData(G,fG(:,4)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));


%% level j=2 basis vectors
BS = LevelBasisSpec(GP,2);
drow = [2,rs(2,3)+2,rs(3,3)+3,rs(4,3)+4];
dvec = zeros(N,4);
dvec(drow(1),1) = 1;
dvec(drow(2),2) = 1;
dvec(drow(3),3) = 1;
dvec(drow(4),4) = 1;
fG =  GHWT_Synthesis(dvec,GP,BS);

[j,k,l] = GHWT_jkl(GP,drow(1),3);
G9 = EditName(ReplaceData(G,fG(:,1)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = GHWT_jkl(GP,drow(2),3);
G10 = EditName(ReplaceData(G,fG(:,2)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = GHWT_jkl(GP,drow(3),3);
G11 = EditName(ReplaceData(G,fG(:,3)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));

[j,k,l] = GHWT_jkl(GP,drow(4),3);
G12 = EditName(ReplaceData(G,fG(:,4)),sprintf('$$\\psi^{%d}_{%d,%d}$$',j,k,l));


%% generate GHWT figures in Matlab

g = figure;

s1 = subplot(3,4,1);
axis equal;
axis off;
cmax = max(abs(G1));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G1);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G1);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,2);
axis equal;
axis off;
cmax = max(abs(G2));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G2);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G2);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,3);
axis equal;
axis off;
cmax = max(abs(G3));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G3);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G3);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,4);
axis equal;
axis off;
cmax = max(abs(G4));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G4);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G4);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,5);
axis equal;
axis off;
cmax = max(abs(G5));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G5);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G5);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,6);
axis equal;
axis off;
cmax = max(abs(G6));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G6);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G6);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,7);
axis equal;
axis off;
cmax = max(abs(G7));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G7);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G7);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,8);
axis equal;
axis off;
cmax = max(abs(G8));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G8);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G8);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,9);
axis equal;
axis off;
cmax = max(abs(G9));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G9);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G9);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,10);
axis equal;
axis off;
cmax = max(abs(G10));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G10);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G10);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,11);
axis equal;
axis off;
cmax = max(abs(G11));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G11);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G11);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

s1 = subplot(3,4,12);
axis equal;
axis off;
cmax = max(abs(G12));
set(gca, 'CLim', [-cmax, cmax]);
[~,~,~,name] = ExtractData(G12);
title(name,'interpreter','latex','FontSize',16);
GraphSig_Plot(G12);
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

PartitionTreeDisplay(GP,1);
title('GHWT fine-to-coarse dictionary structure');

[c2f,f2c] = HGLET_GHWT_BasisCount(GP);
c2f = sprintf('10^%d = 2^%d',floor(-c2f),floor(-c2f/log10(2)));
f2c = sprintf('10^%d = 2^%d',floor(-f2c),floor(-f2c/log10(2)));

fprintf('\n\nThere are (rounded due to machine arithmetic):\n\n');
fprintf('* %s choosable orthonormal bases from the GHWT coarse-to-fine dictionary\n',c2f);
fprintf('* %s choosable orthonormal bases from the GHWT fine-to-coarse dictionary\n                  (not contained in the coarse-to-fine dictionary)\n',f2c);