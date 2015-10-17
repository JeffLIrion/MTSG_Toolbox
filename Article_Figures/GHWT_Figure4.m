% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% NOTE: These figures and results will differ from those in the article
% because the noise is generated randomly and thus not the same as in the
% experiment that produced the article figures.  

load('MN_MutGauss.mat');

% convert to inverse Euclidean weights
G = Adj2InvEuc(G);

% partition the graph
GP = PartitionTreeFiedler(G,1);

% add noise to the mutilated Gaussian
GN = AddNoise(G,5,'gaussian');

% perform the GHWT transform on the noisy signal
dmatrix = GHWT_Analysis(GN,GP);

% find the GHWT best basis
[dvec,BS] = GHWT_BestBasis(dmatrix,GP,0.1);

% threshold the GHWT coefficients
[dvecT,kept] = dvec_Threshold(dvec,'s',0.11,GP,BS);

% reconstruct
[~,G4c] = GHWT_Synthesis(dvecT,GP,BS,GN);

% print the SNR of the reconstructed signal
fprintf('\n\nThe GHWT reconstructed signal has SNR %.2f dB (%d%% coefficients kept)\n\n',SNR(G,G4c),round(100*kept/length(dvecT)));

% take the absolute value of the original minus the denoised signal
G4d = abs(GN-G4c);

% find the min and max values for the colormap
cmin = min([min(G), min(GN), min(G4c), min(G4d)]);
cmax = max([max(G), max(GN), max(G4c), max(G4d)]);

% set the min and max values for the colormap
G4a = EditPlotSpecs(G,sprintf('notitle CLim[%f,%f]',cmin,cmax));
G4b = EditPlotSpecs(GN,sprintf('notitle CLim[%f,%f]',cmin,cmax));
G4c = EditPlotSpecs(G4c,sprintf('notitle CLim[%f,%f]',cmin,cmax));
G4d = EditPlotSpecs(G4d,sprintf('notitle CLim[%f,%f]',cmin,cmax));


%%% display the figures

figure

% Figure 4a
s1 = subplot(2,2,1);
axis equal;
axis off;
set(gca, 'CLim', [cmin, cmax]);
title('Original Signal');
GraphSig_Plot(G4a);

% grab the plot and put it in the figure with all the basis vectors
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% Figure 4b
s1 = subplot(2,2,2);
axis equal;
axis off;
set(gca, 'CLim', [cmin, cmax]);
title('Noisy Signal');
GraphSig_Plot(G4b);

% grab the plot and put it in the figure with all the basis vectors
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% Figure 4c
s1 = subplot(2,2,3);
axis equal;
axis off;
set(gca, 'CLim', [cmin, cmax]);
title('Denoised via GHWT');
GraphSig_Plot(G4c);

% grab the plot and put it in the figure with all the basis vectors
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

% Figure 4d
s1 = subplot(2,2,4);
axis equal;
axis off;
set(gca, 'CLim', [cmin, cmax]);
title('abs(Original-Denoised)');
GraphSig_Plot(G4d);

% grab the plot and put it in the figure with all the basis vectors
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


set(gcf,'color','w');