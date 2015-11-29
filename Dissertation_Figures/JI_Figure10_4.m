% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



matrix = imread('barbara.png','png');
% matrix = imresize(matrix,0.5);
matrix = double(matrix)/255;

% constants
[rows,cols] = size(matrix);

% generate row and column permutation vectors
[~,Prows] = sort(rem((1:rows)'*pi,1),'descend');        % the row permutation
[~,Pcols] = sort(rem((1:cols)'*sqrt(2),1),'descend');   % the reverse row permutation
rowsP = (1:rows)'; rowsP(Prows) = rowsP;                % the reverse row permutation
colsP = (1:cols)'; colsP(Pcols) = colsP;                % the reverse col permutation

% permute the matrix
matrix = matrix(Prows,Pcols);
figure;
imshow(matrix);
colormap(gray);
set(gcf,'color','w');

% partition the matrix
[GProws,GPcols] = PartitionTreeMatrixDhillon(matrix);

% specify global row and column bases ==> no partitions will be drawn
BSrows = LevelBasisSpec(GProws,0);
BScols = LevelBasisSpec(GPcols,0);

% display the recursive partitionings
Matrix_BasisVisual(matrix,GProws,GPcols,BSrows,BScols);