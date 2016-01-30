% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



load('Dendrite.mat','G');

% the number of vertices
N = length(G);

% partition the graph
GP = PartitionTreeFiedler(G);


%% for this figure we use ~half of the graph (cut by the Fiedler vector)

psi = GHWT_BasisVector(G,GP,2,1);
[W,xy,f] = ExtractData(psi);
for node = N:-1:1
    if f(node) > 0
        W(node,:) = [];
        W(:,node) = [];
        f(node) = [];
        xy(node,:) = [];
    end
end
N = length(f);

% rotate the graph by t radians
t = pi/4;
R = [cos(t), -sin(t); sin(t), cos(t)];
xy = xy*R';
G = GraphSig(W,xy,f);
GP = PartitionTreeFiedler(G);
[ind,rs] = ExtractData(GP);
[~,jmax] = size(rs);
jmax=jmax-1;


%% GHWT scaling vectors

psi_j0_k0 = GHWT_BasisVector(G,GP,rs(1,1),1); psi_j0_k0 = psi_j0_k0/norm(psi_j0_k0,inf);

psi_j1_k0 = GHWT_BasisVector(G,GP,rs(1,2),2); psi_j1_k0 = psi_j1_k0/norm(psi_j1_k0,inf);
psi_j1_k1 = GHWT_BasisVector(G,GP,rs(2,2),2); psi_j1_k1 = psi_j1_k1/norm(psi_j1_k1,inf);

psi_j2_k0 = GHWT_BasisVector(G,GP,rs(1,3),3); psi_j2_k0 = psi_j2_k0/norm(psi_j2_k0,inf);
psi_j2_k1 = GHWT_BasisVector(G,GP,rs(2,3),3); psi_j2_k1 = psi_j2_k1/norm(psi_j2_k1,inf);
psi_j2_k2 = GHWT_BasisVector(G,GP,rs(3,3),3); psi_j2_k2 = psi_j2_k2/norm(psi_j2_k2,inf);
psi_j2_k3 = GHWT_BasisVector(G,GP,rs(4,3),3); psi_j2_k3 = psi_j2_k3/norm(psi_j2_k3,inf);

% the weight matrix of the cut graph at level j=1
W1 = sparse(N,N);
indrs = ind(rs(1,2):rs(2,2)-1); W1(indrs,indrs) = W1(indrs,indrs) + W(indrs,indrs);
indrs = ind(rs(2,2):rs(3,2)-1); W1(indrs,indrs) = W1(indrs,indrs) + W(indrs,indrs);

% the weight matrix of the cut graph at level j=2
W2 = sparse(N,N);
indrs = ind(rs(1,3):rs(2,3)-1); W2(indrs,indrs) = W2(indrs,indrs) + W(indrs,indrs);
indrs = ind(rs(2,3):rs(3,3)-1); W2(indrs,indrs) = W2(indrs,indrs) + W(indrs,indrs);
indrs = ind(rs(3,3):rs(4,3)-1); W2(indrs,indrs) = W2(indrs,indrs) + W(indrs,indrs);
indrs = ind(rs(4,3):rs(5,3)-1); W2(indrs,indrs) = W2(indrs,indrs) + W(indrs,indrs);

% generate GraphSig objects for levels j=0, j=1, j=2, and j=jmax
Level_0 = psi_j0_k0;
Level_1 = GraphSig(W1,xy,f); Level_1 = ReplaceData(Level_1,0*psi_j1_k0 + psi_j1_k1);
Level_2 = GraphSig(W2,xy,f); Level_2 = ReplaceData(Level_2,0*psi_j2_k0 + 1/3*psi_j2_k1 + 2/3*psi_j2_k2 + psi_j2_k3);
Level_jmax = 1*GraphSig(sparse(N,N),xy,rand(N,1));

% modify the plotspecs
Level_0 = EditPlotSpecs(Level_0,'notitle size5 CLim[0,1] nocolorbar verbatim{{axis off; axis equal;}}');
Level_1 = EditPlotSpecs(Level_1,'notitle size5 CLim[0,1] nocolorbar verbatim{{axis off; axis equal;}}');
Level_2 = EditPlotSpecs(Level_2,'notitle size5 CLim[0,1] nocolorbar verbatim{{axis off; axis equal}}');
Level_jmax = EditPlotSpecs(Level_jmax,'notitle size5 CLim[0,1] nocolorbar verbatim{{axis off; axis equal;}}');


%% display the number of choosable bases
[c2f,f2c] = HGLET_GHWT_BasisCount(GP);

fprintf('\n\nThere are:\n\n');
fprintf('* %2d choosable orthonormal bases from the HGLET (L) dictionary\n\n',c2f);
fprintf('* %2d choosable bases from the HGLET (Lrw) dictionary\n\n',c2f);
fprintf('* %2d choosable orthonormal bases from the HGLET (Lsym) dictionary\n\n',c2f);
fprintf('* %2d choosable orthonormal bases from the GHWT coarse-to-fine dictionary\n',c2f);
fprintf('* %2d choosable orthonormal bases from the GHWT fine-to-coarse dictionary\n                (not contained in the coarse-to-fine dictionary)\n',f2c);
fprintf('* %2d total choosable orthonormal bases from the GHWT dictionaries\n\n',c2f+f2c);


%% display the four figures
h = figure;

% plot Level 0
s1 = subplot(1,4,1);
axis equal;
axis off;
set(gca, 'CLim', [0, 1]);
title('Level 0');

GraphSig_Plot(Level_0);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% plot Level 1
s1 = subplot(1,4,2);
axis equal;
axis off;
set(gca, 'CLim', [0, 1]);
title('Level 1');

GraphSig_Plot(Level_1);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% plot Level 2
s1 = subplot(1,4,3);
axis equal;
axis off;
set(gca, 'CLim', [0, 1]);
title('Level 2');

GraphSig_Plot(Level_2);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% plot Level jmax
s1 = subplot(1,4,4);
axis equal;
axis off;
set(gca, 'CLim', [0, 1]);
title(sprintf('Level %d (jmax)',jmax));

GraphSig_Plot(Level_jmax);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


colormap('jet');
set(gcf,'color','w');