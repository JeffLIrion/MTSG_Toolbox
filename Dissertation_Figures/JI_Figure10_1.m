% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~exist('ScienceNews.mat','file')
    fprintf('\n\n\nIn order to generate these results, install Python,\nrun WavePath.m, then run this script again.\n\n\n');
    return
end

load('ScienceNews.mat');
[rows,cols] = size(cmatrix);

% display the original matrix
figure;
spy(cmatrix);
h1 = findobj('Type','line');
set(h1(1),'MarkerSize',1);
axis equal;
axis off;
set(gcf,'color','w');

% partition the matrix
[GProws,GPcols] = PartitionTreeMatrixDhillon(cmatrix);

% display the partitioned matrix
SparseMatrix_BasisVisual(cmatrix,GProws,GPcols);
h1 = findobj('Type','line');
set(h1(3),'MarkerSize',1);

% count the number of choosable bases
[rows_c2f,rows_f2c] = HGLET_GHWT_BasisCount(GProws);
[cols_c2f,cols_f2c] = HGLET_GHWT_BasisCount(GPcols);

fprintf('\n\n');
fprintf('There are %.2e choosable bases for the rows.\n',rows_c2f+rows_f2c);
fprintf('There are %.2e choosable bases for the columns.\n',cols_c2f+cols_f2c);
fprintf('\n\n');