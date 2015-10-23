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
%   sortnodes       plot the signal values from smallest to largest in 
%                   magnitude
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
[temp,CLim] = find_between(G.plotspecs,'CLim[',']');
if CLim
    temp = eval(sprintf('[%s]',temp));
    cmin = temp(1);
    cmax = temp(2);
else
    [temp,CLim] = find_between(G.plotspecs,'CLim([','])');
    if CLim
        temp = eval(sprintf('[%s]',temp));
        cmin = temp(1);
        cmax = temp(2);
    else
        cmin = [];
        cmax = [];
    end
end

% how big should the points be?
if G.dim == 1 
    ptsize = 2;
else
    ptsize = 25;
end

[temp,TF] = find_between(G.plotspecs,'size[',']');
if TF
    ptsize = eval(sprintf('[%s]',temp));
else
    [temp,TF] = find_after(G.plotspecs,'size');
    if TF
        ptsize = temp;
    end
end

% how wide should the lines be?
if G.dim == 1
    linewide = 1;
else
    linewide = [];
end

[temp,TF] = find_after(G.plotspecs,'idth');
if TF
    linewide = str2double(temp);
end

% what color should the lines (1D) / graph edges (2D & 3D) be?
if G.dim == 1
    linecolor = 'b';
else
    linecolor = 'k';
end

[temp,TF] = find_after(G.plotspecs,'olor','olorbar');
if TF
    linecolor = temp;
else
    [temp,TF] = find_between(G.plotspecs,'olor[',']');
    if TF
        linecolor = eval(sprintf('[%s]',temp));
    end
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
[verbtext,verbatim] = find_between(G.plotspecs,'verbatim{{','}}');

% sort the nodes?
sortnodes = false;
if ~isempty(strfind(G.plotspecs,'sortnodes'))
    sortnodes = true;
end


end





function [str,TF] = find_between(plotspecs,left,right)
% Given a string 'plotspecs', return str, where <left>str<right>

str = [];
TF = false;

a = strfind(plotspecs,left);
if ~isempty(a) && length(plotspecs) > a(end)+length(left)-1
    plotspecs(1:a(end)+length(left)-1) = [];
    b = strfind(plotspecs,right);
    if ~isempty(b)
        str = plotspecs(1:b(1)-1);
        TF = true;
    end
end
end




function [y,TF] = find_after(plotspecs,left,mismatch)
% Given a string 'plotspecs', return the number OR character after <left>,
% but don't match 'left' with 'mismatch'

y = [];
TF = false;

a = strfind(plotspecs,left);

if ~isempty(a) && exist('mismatch','var')
    a2 = strfind(plotspecs,mismatch);
    a = setdiff(a,a2);
end

if ~isempty(a) && length(plotspecs) > a(end)+length(left)-1
    if isstrprop(plotspecs(a(end)+length(left)),'alpha')
        y = plotspecs(a(end)+length(left));
        TF = true;
    elseif isstrprop(plotspecs(a(end)+length(left)),'digit')
        onedot = false;
        b = a(end)+length(left);
        while b < length(plotspecs)
            if isstrprop(plotspecs(b+1),'digit')
                b = b+1;
            elseif strcmp(plotspecs(b+1),'.') && ~onedot
                b = b+1;
                onedot = true;
            else
                break
            end
        end
        y = str2double(plotspecs(a(end)+length(left):b));
        TF = true;
    end
end
end