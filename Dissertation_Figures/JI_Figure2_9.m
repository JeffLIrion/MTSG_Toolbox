% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% the grid is Nx by Ny, where K Ny < Nx < (K+1) Ny
K = 10;
Ny = 10;
Nx = K*Ny+1;
N = Nx*Ny;

% construct the graph
G = Ggrid(Nx,Ny);
G = EditPlotSpecs(G,'notitle symm nocolorbar verbatim{{axis equal; axis off;}}');
GP = PartitionTreeFiedler(G);

% specify that we are using the global basis
BS = LevelBasisSpec(GP,0);

% specify the expansion coefficients which will generate the eigenvectors
dvec = zeros(N,2);
dvec(2,1) = 1;
dvec(K+1,2) = 1;
dvec(K+2,3) = 1;

% compute the eigenvectors via HGLET (all at once)
f = HGLET_Synthesis(dvec,GP,BS,G);

% first eigenvector
G1 = ReplaceData(G,f(:,1));

% last x eigenvector
G2 = ReplaceData(G,f(:,2));

% first y eigenvector
G3 = ReplaceData(G,f(:,3));


% display the three figures
h = figure;

% plot G1
s1 = subplot(3,1,1);
axis equal;
axis off;
cmax = max(abs(G1));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\phi_{1}$$'),'interpreter','latex','FontSize',16);

GraphSig_Plot(G1);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% plot G2
s1 = subplot(3,1,2);
axis equal;
axis off;
cmax = max(abs(G2));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\phi_{%d}$$',K),'interpreter','latex','FontSize',16);

GraphSig_Plot(G2);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% plot G3
s1 = subplot(3,1,3);
axis equal;
axis off;
cmax = max(abs(G3));
set(gca, 'CLim', [-cmax, cmax]);
title(sprintf('$$\\phi_{%d}$$',K+1),'interpreter','latex','FontSize',16);

GraphSig_Plot(G3);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


colormap('jet');
set(gcf,'color','w');