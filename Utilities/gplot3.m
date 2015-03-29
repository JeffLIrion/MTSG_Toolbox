function [X,Y,Z] = gplot3(W,xyz,lc,varargin)
% A modification of gplot which can handle 2-D or 3-D graphs and allows
% more customization of the resulting plot.  
%
% Input
%   W           the weight matrix
%   xyz         the matrix of xyz coordinates
%   lc          a line specification that determines line style, marker
%               symbol, and color of the plotted lines
%   varargin    specifications such as the line width
%
% Output
%   X           x-coordinates of the line segments
%   Y           y-coordinates of the line segments
%   Z           z-coordinates of the line segments
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% find the endpoints of the edges ==> the edges are (i,j)
[i,j] = find(W);
[~, p] = sort(max(i,j));
i = i(p);
j = j(p);

% find out whether the graph is 2-D or 3-D
[~,cols] = size(xyz);
if cols == 1
    xyz = [xyz, zeros(length(xyz),1)];
    cols = 2;
elseif cols ~= 2 && cols ~= 3
    return
end

% Create a list of the line segments
X = [ xyz(i,1) xyz(j,1) NaN(size(i))]';
Y = [ xyz(i,2) xyz(j,2) NaN(size(i))]';
X = X(:);
Y = Y(:);
if cols == 3
    Z = [ xyz(i,3) xyz(j,3) NaN(size(i))]';
    Z = Z(:);
end

if nargout == 0
    % specify 'lc" if it is not given as an input
    if ~exist('lc','var') || ~ischar(lc)
        lc = '-k';
    end
    
    % plot the graph
    if cols ==2
        h = plot(X,Y,lc);
    elseif cols == 3
        h = plot3(X(:),Y(:),Z(:),lc);
    end
    
    % customize the plot
    for k=1:nargin-3
        s = varargin{k};
        try
            set(h,s,varargin{k+1});
        catch

        end
    end
end
