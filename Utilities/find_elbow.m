function IX = find_elbow(x,y)
% Given the x and y coordinates of a curve, find the elbow.  
% Inspired by http://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
%
% Input
%   x       the x coordinates of the points
%   y       the y coordinates of the points
%
% Output
%   IX      the index of the elbow
%
%
% 
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% a unit vector pointing from (x1,y1) to (xN,yN)
v = [x(end)-x(1), y(end)-y(1)];
v = v/norm(v,2);

% subtract (x1,y1) from the coordinates
xy = [x-x(1), y-y(1)];

% the hypotenuse
H = (sum(xy.^2,2)).^0.5;

% the adjacent side
A = xy*v';

% the opposite side
O = (H.^2-A.^2).^0.5;

% find the largest distance
[~,IX] = max(O);


end