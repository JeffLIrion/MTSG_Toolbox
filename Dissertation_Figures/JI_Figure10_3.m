% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



matrix = imread('barbara.png','png');
% matrix = imresize(matrix,0.5);
matrix = double(matrix);

figure;
imshow(matrix/255);
colormap(gray);
set(gcf,'color','w');