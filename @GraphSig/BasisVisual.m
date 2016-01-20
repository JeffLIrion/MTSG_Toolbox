function [dtable,G_regions] = BasisVisual(G,GP,BS,dvec,edge_color)
% Display an HGLET or GHWT basis
%
% Input
%   G           a GraphSig object
%   GP          a GraphPart object
%   BS          a BasisSpec object
%   dvec        the expansion coefficients corresponding to the chosen
%               basis
%   edge_color  the color for displaying the edges that are cut for an
%               HGLET or GHWT coarse-to-fine basis
%
% Output
%   dtable      a table of the specified basis coefficients (use 'imagesc')
%   G_regions   a GraphSig object that displays the regions of the 
%               specified basis
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% extract data and figure out which dictionary I'm dealing with
[ind,rs] = ExtractData(GP);
[N,jmax] = size(rs);
N = N-1;
[levlist,levlengths,BSc2f] = ExtractData(BS,GP);


%% 1. The table of specified coefficients

if exist('dvec','var')
    dtable = -max(abs(dvec))*ones(jmax,N)/63;
else
    dtable = zeros(jmax,N);
end

indr0 = 1;
for row = 1:length(levlist)
    indr = indr0:indr0+levlengths(row)-1;
    if exist('dvec','var') == 1
        dtable(levlist(row),indr) = abs(dvec(indr));
    else
        dtable(levlist(row),indr) = 1;
    end
    indr0 = indr0+levlengths(row);
end

% display the table of specified coefficients
figure;
imagesc(dtable);
if exist('dvec','var')
    colormap('jet');
    cmap = colormap;
    cmap(1,:) = 1;
    colormap(cmap);
    colorbar;
else
    cmap = ones(64,3);
    cmap(end,:) = [0, 0, 0.5];
    colormap(cmap);
end

% labeling & formatting
xlabel('Coefficient Index');
ylabel('Level (j)');
ylim([0.5,jmax+0.5]);

% the tick marks
YT = get(gca,'YTick');
YT = unique([0,round(YT)]);
YT(YT < 1 | YT > jmax) = [];
set(gca,'YTick',YT);

% the tick mark labels
YTL = get(gca,'YTickLabel');
if iscell(YTL)
    YTL = str2num(char(YTL));
else
    YTL = str2num(YTL);
end
if BSc2f
    title('Coefficients of the Specified Basis');
    set(gca,'YTickLabel',YTL-1);
else
    set(gca,'YTickLabel',jmax-YTL);
    title('Coefficients of the Specified Basis (Fine-To-Coarse Dictionary)');
end

set(gcf,'color','w');


%% 2. The GraphSig object for GraphSig_Overlay (coarse-to-fine only)

if BSc2f && G.dim > 0 && G.dim < 4
    G_regions = G;
    G_regions.W = sparse(N,N);

    % only retain edge weights kept intact in the basis' regions
    indr0 = 1;
    for row = 1:length(levlist)
        indr = indr0:indr0+levlengths(row)-1;
        inds = ind(indr);
        G_regions.W(inds,inds) = G.W(inds,inds);
        indr0 = indr0+levlengths(row);
    end

    % the signal on the GraphSig object ==> region levels
    G_regions.f = zeros(N,1);
    G_regions.f(ind) = BSfull(GP,BS)-1;

    % display the overlaid GraphSig object
    G_regions = EditName(G_regions,'Regions of the Specified Basis');
    G_regions = EditPlotSpecs(G_regions,sprintf(' CLim[%d,%d]',0,jmax-1),1);
    if ~exist('edge_color','var')
        edge_color = [255, 102, 178]/255;
    end
    GraphSig_Overlay(G_regions,G,1,'-','k',1,'-',edge_color);
    
    % the colorbar tick marks
    YTC = get(colorbar,'YTick');
    YTC = unique(round([0,YTC]));
    YTC( YTC < 0 | YTC > jmax-1) = [];
    set(colorbar,'YTick',YTC);
else
    G_regions = [];
end


end