% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



fprintf('\n\nThis script does not generate the approximation results for the Graph-QMF and the Laplacian multiwavelets.\n\n');


%% 0. Preliminaries
load('Toronto.mat','G');
N = length(G);

% partition the graph
GP = PartitionTreeFiedler(G,1);


%% 1. Analyze the signal

[dmatrixH,dmatrixHrw,dmatrixHsym] = HGLET_Analysis_All(G,GP);
dmatrixG = GHWT_Analysis(G,GP);

% 1) HGLET best basis with L
[dvec1,BS1,~,p1] = HGLET_GHWT_BestBasis_minrelerror(dmatrixH,0,0,0,GP,G);
r1 = orth2relerror(dvec1);

% 2) Laplacian eigenvectors
BS2 = LevelBasisSpec(GP,0);
dvec2 = dmatrix2dvec(dmatrixH,GP,BS2);
r2 = orth2relerror(dvec2);

% 3) GHWT best basis
[dvec3,BS3,~,p3] = HGLET_GHWT_BestBasis_minrelerror(0,0,0,dmatrixG,GP,G);
r3 = orth2relerror(dvec3);

% 4) Haar basis
BS4 = HaarBasisSpec(GP);
dvec4 = dmatrix2dvec(dmatrixG,GP,BS4);
r4 = orth2relerror(dvec4);

% 5) Hybrid best basis (HGLET with L/Lrw/Lsym + GHWT)
[dvec5,BS5,trans5,p5] = HGLET_GHWT_BestBasis_minrelerror(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GP,G,1);
B5 = HGLET_GHWT_Synthesis(eye(N),GP,BS5,trans5,G);
r5 = nonorth2relerror(dvec5,B5);

% 6) Walsh basis
BS6 = LevelBasisSpec(GP,0);
dvec6 = dmatrix2dvec(dmatrixG,GP,BS6);
r6 = orth2relerror(dvec6);


%% Figure 4
x = (1:N)'/N;

% add eps to the relative errors to avoid zeros
r1=r1+eps;
r2=r2+eps;
r3=r3+eps;
r4=r4+eps;
r5=r5+eps;
r6=r6+eps;

% for generating the legend cell
legcell = {};
row = 0;

figure;
semilogy(x,r6,':m','LineWidth',3); row=row+1; legcell{row} = 'Walsh';
hold on
semilogy(x,r2,'-b','LineWidth',3); row=row+1; legcell{row} = 'Laplacian Eigenvectors (L)';
semilogy(x,r4,'-g','LineWidth',3); row=row+1; legcell{row} = 'Haar';
semilogy(x,r1,':','Color',[0 0.5 0],'LineWidth',3); row=row+1; legcell{row} = sprintf('HGLET (L) BB (tau=%.1f)',p1);
semilogy(x,r3,'-','Color',[255 129 0]/255,'LineWidth',3); row=row+1; legcell{row} = sprintf('GHWT BB (tau=%.1f)',p3);
semilogy(x,r5,'-','Color',[0.8 0.8 0],'LineWidth',3); row=row+1; legcell{row} = sprintf('Hybrid BB (tau=%.1f)',p5);

xlabel('Fraction of Coefficients Retained');
ylabel('Relative Approximation Error');
legend(legcell,'Location','southwest');

xlim([0,1]);
ylim([10^-8,1]);
set(gcf,'color','w');