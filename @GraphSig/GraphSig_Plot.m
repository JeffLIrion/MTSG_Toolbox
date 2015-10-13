function fig = GraphSig_Plot(G)
% Display a plot of the data in a GraphSig object
%
% Input
%   G       a GraphSig object
%
% Output
%    fig    a plot figure
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% only run this for 1D signals
[~,fcols] = size(G.f);
if fcols == 0
    G.f = zeros(G.length,1);
elseif fcols ~= 1
    return
end


%% Case 1: the graph does not have spatial coordinates or it is too large
if G.length > 10^4 || G.dim > 3 || G.dim < 1
    if G.length < 10^6
        fig = figure('visible','on');
        spy(G.W,'o');
        markerH = findall(gca,'color','b');
        set(markerH,'Color','k');
        set(markerH,'MarkerFaceColor','k');
        ezis = get(markerH,'MarkerSize');
        hold on
        scatter((1:G.length)',(1:G.length)',4*ezis,G.f,'fill');
        title(G.name);
    else
        fig = 0;
    end
    return
end


%% Case 2: the graph does have spatial coordinates

% extract the plot specifications
[symmetric,graybar,gray255bar,copperbar,notitle,nocolorbar,stemplot,CLim,cmin,cmax,ptsize,linewide,linecolor,marker,verbatim,verbtext,sortnodes] = ExtractPlotSpecs(G);
set(0, 'DefaultFigureVisible', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 1-D %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if G.dim == 1
    fig = figure('visible','on');
    if G.length < 65 || stemplot
        stem(G.xy, G.f, '-o','LineWidth',2,'Color',linecolor,'MarkerFaceColor',linecolor,'MarkerSize',ptsize(1));
    else
        plot(G.xy, G.f, '-','Color',linecolor,'LineWidth',linewide);
    end

    xmin = min(G.xy);
    xmax = max(G.xy);
    
%     dx = xmax-xmin;
%     if dx < 10*eps
%         dx = 1;
%     end
%     xmin = xmin - 0.1*dx;
%     xmax = xmax + 0.1*dx;

    ymin = min(G.f);
    ymax = max(G.f);
    dy = ymax-ymin;
    if dy < 10*eps
        dy = 1;
    end
    ymin = ymin - 0.1*dy;
    ymax = ymax + 0.1*dy;

    % symmetric?
    if symmetric
        ymax = 1.1*max(ymax, -ymin);
        ymin = -ymax;
    end

    axis([xmin, xmax, ymin, ymax]);
    
    if CLim
        axis([xmin, xmax, cmin, cmax]);
    else
        axis([xmin, xmax, ymin, ymax]);
    end

    % title?
    if ~isempty(G.name) && ~notitle
        title(G.name);
    end
    
    % verbatim?
    if verbatim
        eval(verbtext);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 2-D  or 3-D %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    
    fig = figure('visible','on');
    
    if sortnodes
        [~,IX] = sort(abs(G.f),'descend');
        G.W = G.W(IX,IX);
        G.xy = G.xy(IX,:);
        G.f = G.f(IX);
    end
    
    % plot the graph
    if isempty(linewide)
        gplot3(G.W,G.xy,'LineStyle','-','Color',linecolor);
    else
        gplot3(G.W,G.xy,'LineStyle','-','Color',linecolor,'LineWidth',linewide);
    end
    hold on
    
    % plot the nodes
    if ptsize(1) ~= 0
        % if we're using various point sizes, determine the vector of sizes
        if ~isscalar(ptsize)
            if CLim
                ptsize = ptsize(1) + ptsize(2)*abs(G.f)/max(abs([cmin,cmax]));
            else
                ptsize = ptsize(1) + ptsize(2)*abs(G.f)/max(abs(G.f));
            end
        end
        
        % plot the nodes for a 2-D graph
        if G.dim == 2
            if ischar(marker)
                scatter(G.xy(:,1), G.xy(:,2), ptsize, marker, 'filled');
            else
                scatter(G.xy(:,1), G.xy(:,2), ptsize, G.f, 'filled');
            end

        % plot the nodes for a 3-D graph
        else
            if ischar(marker)
                scatter3(G.xy(:,1), G.xy(:,2), G.xy(:,3), ptsize, marker, 'filled');
            else
                scatter3(G.xy(:,1), G.xy(:,2), G.xy(:,3), ptsize, G.f, 'filled');
            end
        end
    
        % colorbar?
        if ~nocolorbar
            colorbar;
        end

        % custom dynamic display range?
        if CLim
            set(gca, 'CLim', [cmin, cmax]);
        end

        % symmetric?
        if symmetric
            if CLim
                cmax = max([abs(cmin), abs(cmax)]);
            else
                cmax = max(abs(G.f));
            end
            if cmax < 10*eps
                cmax = 1;
            end
            set(gca, 'CLim', [-cmax, cmax]);
        end

        % gray255?
        if gray255bar
            colormap('gray');
            set(gca, 'CLim', [0, 255]);

        % gray?
        elseif graybar
            colormap('gray');

        % copper?
        elseif copperbar
            colormap('copper');
        end
    end

    % title?
    if ~isempty(G.name) && ~notitle
        title(G.name);
    end
    
    % verbatim?
    if verbatim
        eval(verbtext);
    end
    
end

set(gcf,'color','w');

set(gca, 'FontSize',10);

% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0.01 0.01 4.77 3.91]);
% set(gcf, 'PaperSize', [4.78 3.92]);
% set(gca, 'FontSize',10);

% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [5 5]);
% set(gcf, 'PaperPosition', [0 0 5 5]);


% if G.dim == 2
%     annotation(fig,'rectangle',[0.08 0.05 0.83 0.89],'LineWidth',0.1,'FaceColor','flat','EdgeColor','w');
% end



end