% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% load the true signal
load('Dendrite.mat');
G0 = G;

% load the noisy signal
load('Dendrite_SNR8.mat');

% partition the graph
GP = PartitionTreeFiedler(G,1);


% analyze the true and noisy signals
dG = GHWT_Analysis(G,GP);
dG0 = GHWT_Analysis(G0,GP);

% find the best basis for the noisy signal
[dvec,BS,~,tau] = HGLET_GHWT_BestBasis_minrelerror(0,0,0,dG,GP,G);
dvec0 = dmatrix2dvec(dG0,GP,BS);


% the magnitudes of the coefficients
x = sort(abs(dvec),'descend');

% compute the relative error curve
r = orth2relerror(dvec);

% compute the SNR curve
s = orth2SNR(dvec,dvec0,GP,BS,'s');

% shift the data points
x = [x; 0];
r = [r(1); r];
s = [s(1); s];

% find the threshold
idx = find_elbow(x,r);
idx = idx-1+find_elbow(x(idx:end),r(idx:end));

% the denoised signal
dvec = dvec_Threshold(dvec,'s',idx,GP,BS);
[~,G2] = GHWT_Synthesis(dvec,GP,BS,G);

% find the color limits
cmax = max([max(G0); max(G); max(G2)]);
cmin = min([min(G0); min(G); min(G2)]);

% prepare to display the true, noisy, and denoised signals
G0 = EditPlotSpecs(G0,sprintf(' notitle CLim([%f,%f])',cmin,cmax),1);
G  = EditPlotSpecs(G, sprintf(' notitle CLim([%f,%f])',cmin,cmax),1);
G2 = EditPlotSpecs(G2,sprintf(' notitle CLim([%f,%f])',cmin,cmax),1);

% display the figures
h = figure;
colormap('jet');

% original
s1 = subplot(2,2,1);
axis equal;
axis off;
set(gca, 'CLim', [cmin, cmax]);
title('Original');

GraphSig_Plot(G0);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

% noisy
s1 = subplot(2,2,2);
axis equal;
axis off;
set(gca, 'CLim', [cmin, cmax]);
title('Noisy');

GraphSig_Plot(G);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);

% denoised
s1 = subplot(2,2,3);
axis equal;
axis off;
set(gca, 'CLim', [cmin, cmax]);
title('Denoised');

GraphSig_Plot(G2);
ax1 = gca;
h1 = gcf;
fig1 = get(ax1,'children');
copyobj(fig1,s1);
close(h1);


% plot the relative error and SNR curves
s1 = subplot(2,2,4);
[hAx,h1,h2] = plotyy(x,r,x,s);
set(h1,'LineWidth',2,'Color','b')
set(h2,'LineWidth',2,'Color',[0, 0.5, 0])
set(hAx(1),'ycolor','b');
set(hAx(2),'ycolor',[0, 0.5, 0]);
xlabel('Threshold','FontSize',12);
ylabel(hAx(1),'Relative Error','FontSize',12);
ylabel(hAx(2),'SNR','FontSize',12);
set(hAx(1),'XLim',[0 x(1)]);
set(hAx(2),'XLim',[0 x(1)]);
set(hAx(1),'YLim',[0 1]);
axis(hAx,'square');
set(hAx(1),'YTick',[0, 0.25, 0.5, 0.75, 1]);
set(gcf,'color','w');

% draw the vertical line
hold(hAx(1),'on');
line([x(idx),x(idx)],[0,1],'Color','r');


set(gcf,'color','w');