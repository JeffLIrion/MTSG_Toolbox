function G = GraphSig(W,xy,f,dataname,plotspecs)
% Constructor for a GraphSig object, a struct containing the fields:
%
%           - name: name of the GraphSig object
%           - plotspecs: special instructions for plotting the object
%           - length: length of the data sequence
%           - size: size of all the data specified
%           - dim: the dimension of the problem (e.g. 1 for 1D, 2 for 2D)
%           - W: edge weight matrix
%           - xy: spatial coordinates (if they are specified)
%           - f: the data sequence
%
% Input
%   W               the edge weight matrix
%   xy              the coordinates of the vertices (not necessarily 2D)
%   f               the values of the function/data at the vertices
%   dataname        a title for the data
%   plotspecs       specifications for plotting:
%                   - symm: use a symmetric colorscale
%                   - gray: use grayscale
%                   - gray255: plot a grayscale image with color bounds
%                     [0,255]
%                   - copper: use a copper color scale
%                   - notitle: don't display the title
%                   - nocolorbar: don't display a colorbar
%                   - stem: use a stem plot
%                   - CLim[cmin,cmax]: set the dynamic display range to
%                     [cmin,cmax]
%                   - size25: set the size of the nodes to 25 (or whatever
%                     value is specified)
%                   - LineWidth2: set the width of the lines to 2 (or 
%                     whatever value is specified)
%                   - LineColorb: set the color of the lines/graph edges to
%                     color 'b' (or whatever color is specified)
%                   - red/black/blue: make the nodes red, black, or blue
%                   - verbatim{{...}}: execute the command in the brackets
%
% Output
%   G               a struct with eight fields: name, plotspecs, length, 
%                   size, dim, W, xy, and f
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% W
if exist('W','var')
    datastruct.W = W;
else
    datastruct.W = [];
end


%% xy
if exist('xy','var') && length(xy) == length(W)
    datastruct.xy = xy;
else
    datastruct.xy = [];
end


%% f
if exist('f','var')
    [frows,fcols] = size(f);
    if frows == length(W)
        datastruct.f = f;
    elseif fcols == length(W)
        datastruct.f = f';
    else
        datastruct.f = [];
    end
else
    datastruct.f = [];
end


%% name
if exist('dataname','var') && ischar(dataname) == 1
    datastruct.name = dataname;
else
    datastruct.name = [];
end


%% plotspecs
if exist('plotspecs','var') && ischar(plotspecs) == 1
    datastruct.plotspecs = plotspecs;
else
    datastruct.plotspecs = [];
end


%% figure out length, dim, and size
datastruct.length = length(datastruct.W);
[~,datastruct.dim] = size(datastruct.xy);
[~,fcols] = size(datastruct.f);
datastruct.size = size(datastruct.W) + [0,datastruct.dim] + [0,fcols];



G = class(datastruct,'GraphSig');


end