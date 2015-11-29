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

% partition the matrix
[GProws,GPcols] = PartitionTreeMatrixDhillon(matrix);

% "normal" best basis
[dvec,BSrows,BScols] = GHWT_Matrix_Analysis_BestBasis(matrix,GProws,GPcols,1,1,1,1);
relerror = orth2relerror(dvec(:));

% Haar basis
dvec_Haar = GHWT_Matrix_Analysis(matrix,GProws,GPcols,HaarBasisSpec(GProws),HaarBasisSpec(GPcols));
relerror_Haar = orth2relerror(dvec_Haar(:));

% Walsh basis
dvec_Walsh = GHWT_Matrix_Analysis(matrix,GProws,GPcols,LevelBasisSpec(GProws,0),LevelBasisSpec(GPcols,0));
relerror_Walsh = orth2relerror(dvec_Walsh(:));

if exist('WavePath.m','file')
    % extract the indices from the graph partitionings
    ind_rows = ExtractData(GProws);
    ind_cols = ExtractData(GPcols);
    
    % Haar wavelet transform parameters
    wavelet_type = 'Haar';
    coarse_level = 0;
    qmf = MakeONFilter(wavelet_type);

    % 1) Haar wavelet transform on shuffled image
    dvec_Haar1 = FWT2_PO(matrix,coarse_level,qmf);
    relerror_Haar1 = orth2relerror(dvec_Haar1(:));

    % 2) Haar wavelet transform on GP-ordered matrix
    dvec_Haar2 = FWT2_PO(matrix(ind_rows,ind_cols),coarse_level,qmf);
    relerror_Haar2 = orth2relerror(dvec_Haar2(:));
    
    % Coiflet wavelet transform parameters
    wavelet_type = 'Coiflet';
    parameter = 4;
    coarse_level = 0;
    qmf = MakeONFilter(wavelet_type,parameter);

    % 1) Coiflet wavelet transform on shuffled image
    dvec_Coiflet1 = FWT2_PO(matrix,coarse_level,qmf);
    relerror_Coiflet1 = orth2relerror(dvec_Coiflet1(:));

    % 2) Coiflet wavelet transform on GP-ordered matrix
    dvec_Coiflet2 = FWT2_PO(matrix(ind_rows,ind_cols),coarse_level,qmf);
    relerror_Coiflet2 = orth2relerror(dvec_Coiflet2(:));
end

% plot the relative errors
close all
x = 1:length(relerror);

if exist('WavePath.m','file')
    figure;
    plot(x,relerror_Coiflet1,'--r','LineWidth',2);
    hold on
    plot(x,relerror_Coiflet2,'-r','LineWidth',2);
    plot(x,relerror_Haar1,'--','Color',[0, 0.5, 0],'LineWidth',2);
    plot(x,relerror_Haar2,'-','Color',[0, 0.5, 0],'LineWidth',2);
    plot(x,relerror_Haar,'-g','LineWidth',3);
    plot(x,relerror_Walsh,'-c','LineWidth',2);
    plot(x,relerror,':b','LineWidth',2);
    legend('Coiflet-4 (Shuffled)','Coiflet-4 (Reordered)','Haar Wavelet (Shuffled)','Haar Wavelet (Reordered)','Graph Haar','Graph Walsh','GHWT BB (fine-to-coarse)');
else
    fprintf('\n\nPlease install WaveLab if you would like to generate Haar wavelet and Coiflet relative error curves.\n\n<a href="http://statweb.stanford.edu/~wavelab/">http://statweb.stanford.edu/~wavelab/</a>\n\n');
    figure;
    plot(x,relerror_Haar,'-g','LineWidth',3);
    hold on
    plot(x,relerror_Walsh,'-c','LineWidth',2);
    plot(x,relerror,':b','LineWidth',2);
    legend('Haar','Graph Walsh','GHWT BB (fine-to-coarse)');
end
xlabel('Number of Coefficients Retained');
ylabel('Relative Error');
set(gcf,'color','w');