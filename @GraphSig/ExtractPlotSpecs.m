function [symmetric,graybar,gray255bar,copperbar,notitle,nocolorbar,stemplot,CLim,cmin,cmax,ptsize,linewide,linecolor,marker,verbatim,verbtext,sortnodes] = ExtractPlotSpecs(G)
% Extract the important information contained in the 'plotspecs' of a 
% GraphSig object
%
% Input
%   G               a GraphSig object
%
% Output
%   symmetric       symmetrize the colorbar
%   graybar         plot a grayscale image
%   gray255bar      plot a grayscale [0,255] image
%   copperbar       plot a copper image
%   notitle         display a title
%   nocolorbar      display a colorbar
%   stemplot        use a stem plot
%   CLim            specify the dynamic display range
%   cmin            the min of the dynamic display range
%   cmax            the max of the dynamic display range
%   ptsize          the size of the nodes
%   linewide        the width of the lines in gplot
%   linecolor       the color of the lines (1D) / graph edges (2D & 3D)
%   marker          the color of the line (1D) / nodes (2D & 3D)
%   verbatim        use other special instructions
%   verbtext        the specified special instructions
%   
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% symmetrize the colorbar?
symmetric = false;
if ~isempty(strfind(G.plotspecs,'symm'))
    symmetric = true;
end

% plot a grayscale image?
graybar = false;
if ~isempty(strfind(G.plotspecs,'gray'))
    graybar = true;
end

% plot a grayscale [0,255] image?
gray255bar = false;
if ~isempty(strfind(G.plotspecs,'gray255'))
    gray255bar = true;
end

% plot a copper image?
copperbar = false;
if ~isempty(strfind(G.plotspecs,'copper'))
    copperbar = true;
end

% display a title?
notitle = false;
if ~isempty(strfind(G.plotspecs,'notitle'))
    notitle = true;
end

% display a colorbar?
nocolorbar = false;
if ~isempty(strfind(G.plotspecs,'nocolorbar'))
    nocolorbar = true;
end

% use a stem plot?
stemplot = false;
if ~isempty(strfind(G.plotspecs,'stem'))
    stemplot = true;
end

% what should the dynamic display range be?
CLim = false;
s = strfind(G.plotspecs,'CLim[');
if ~isempty(s)
    [cmin,cmax] = bracket_comma_bracket(G.plotspecs(s(end):end));
    if ~isempty(cmin) && ~isempty(cmax)
        CLim = true;
    end
else
    cmin = [];
    cmax = [];
end

% how big should the points be?
if G.dim == 1 
    ptsize = 2;
else
    ptsize = 25;
end
s = strfind(G.plotspecs,'size[');
if ~isempty(s)
    [x1,x2] = bracket_comma_bracket(G.plotspecs(s(end):end));
    if ~isempty(x1) && ~isempty(x2)
        ptsize = [x1;x2];
    end
else
    s = strfind(G.plotspecs,'size');
    l = length(G.plotspecs);
    if ~isempty(s) && l > s(end)+3 && isstrprop(G.plotspecs(s(end)+4),'digit') == 1
        t = s(end)+4;
        onedot = false;
        while t < l
            if isstrprop(G.plotspecs(t+1),'digit')
                t = t+1;
            elseif strcmp(G.plotspecs(t+1),'.') && ~onedot
                t = t+1;
                onedot = true;
            else
                break
            end            
        end
        ptsize = str2double(G.plotspecs(s(end)+4:t));
    end
end

% how wide should the lines be?
l = length(G.plotspecs);
s = strfind(G.plotspecs,'idth');
if ~isempty(s) && l > s(end)+3 && isstrprop(G.plotspecs(s(end)+4),'digit') == 1
    t = s(end)+4;
    onedot = false;
    while t < l
        if isstrprop(G.plotspecs(t+1),'digit')
            t = t+1;
        elseif strcmp(G.plotspecs(t+1),'.') && ~onedot
            t = t+1;
            onedot = true;
        else
            break
        end
    end
    linewide = str2double(G.plotspecs(s(end)+4:t));
else
    if G.dim == 1
        linewide = 1;
    else
        linewide = [];
    end
end

% what color should the lines (1D) / graph edges (2D & 3D) be?
l = length(G.plotspecs);
s = strfind(G.plotspecs,'olor'); %olorbar
while ~isempty(s) && s(end) > 3 && s(end)+5 < l && strcmp(G.plotspecs(s(end)-3:s(end)+6),'nocolorbar')
    s(end) = [];
end
if ~isempty(s) && l > s(end)+3 && isstrprop(G.plotspecs(s(end)+4),'alpha') == 1
    linecolor = G.plotspecs(s(end)+4);
elseif G.dim == 1
    linecolor = 'b';
else
    linecolor = 'k';
end

% what color should the nodes be?
if G.dim == 1
    if ~isempty(strfind(G.plotspecs,'red'))
        marker = '-r';
    elseif ~isempty(strfind(G.plotspecs,'blue'))
        marker = '-b';
    elseif ~isempty(strfind(G.plotspecs,'black'))
        marker = '-k';
    else
        marker = '-b';
    end
    if stemplot
        marker = strcat(marker,'o');
    end
else
    if ~isempty(strfind(G.plotspecs,'red'))
        marker = 'ro';
    elseif ~isempty(strfind(G.plotspecs,'blue'))
        marker = 'bo';
    elseif ~isempty(strfind(G.plotspecs,'black'))
        marker = 'ko';
    else
        marker = [];
    end
end

% other plot instructions?
verbatim = false;
verbstart = strfind(G.plotspecs,'verbatim{{');
if ~isempty(verbstart)
    verbstart = verbstart(1);
    verbend = strfind(G.plotspecs(verbstart+10:end),'}}');
    if ~isempty(verbend)
        verbend = verbend(1)+verbstart+9;
        verbtext = G.plotspecs(verbstart+10:verbend-1);
        verbatim = true;
    end
else
    verbtext = [];
end

% sort the nodes?
sortnodes = false;
if ~isempty(strfind(G.plotspecs,'sortnodes'))
    sortnodes = true;
end


end





function [x1,x2] = bracket_comma_bracket(s)
% Given a string s, return numbers x1 and x2 from s = '...[x1,x2]...'

% the index of '['
bracket1 = strfind(s,'[');
bracket1 = bracket1(1);

% the index of ','
comma = strfind(s,',');
comma = comma(1);

% the index of ']'
bracket2 = strfind(s,']');
bracket2 = bracket2(1);

if bracket1 < comma && comma < bracket2
    try
        x1 = str2double(s(bracket1+1:comma-1));
        x2 = str2double(s(comma+1:bracket2-1));
    catch
        x1 = [];
        x2 = [];
    end
else
    x1 = [];
    x2 = [];
end
end