% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~exist('ScienceNews.mat','file')
    fprintf('\n\n\nIn order to generate these results, install Python,\nrun WavePath.m, then run this script again.\n\n\n');
    return
end

load('ScienceNews.mat');
[rows,cols] = size(cmatrix);

% partition the matrix
[GProws,GPcols] = PartitionTreeMatrixDhillon(cmatrix);
ind_rows = ExtractData(GProws);
ind_cols = ExtractData(GPcols);

% find the GHWT best basis
[dvec,BSrows,BScols] = GHWT_Matrix_Analysis_BestBasis(cmatrix,GProws,GPcols,1,1,1,1);

% the Haar basis
BSrows_Haar = HaarBasisSpec(GProws);
BScols_Haar = HaarBasisSpec(GPcols);
dvec_Haar = GHWT_Matrix_Analysis(cmatrix,GProws,GPcols,BSrows_Haar,BScols_Haar);

% the Walsh basis
BSrows_Walsh = LevelBasisSpec(GProws,0);
BScols_Walsh = LevelBasisSpec(GPcols,0);
dvec_Walsh = GHWT_Matrix_Analysis(cmatrix,GProws,GPcols,BSrows_Walsh,BScols_Walsh);

% generate relative error curves
relerror_GHWT = orth2relerror(dvec(:));
relerror_Haar = orth2relerror(dvec_Haar(:));
relerror_Walsh = orth2relerror(dvec_Walsh(:));

% plot the relative error curves
x = (1:rows*cols)';
figure;
plot(x,relerror_Haar,'-k','LineWidth',2);
hold on
plot(x,relerror_Walsh,'-b','LineWidth',2);
plot(x,relerror_GHWT,'-','Color',[0 0.5 0],'LineWidth',2);
legend('Haar','Walsh','GHWT');
plot([nnz(cmatrix),nnz(cmatrix)],[0,1],'--r','LineWidth',1);
xlim([0,rows*cols+1]);
xlabel('Number of Coefficients Retained');
ylabel('Relative Error');
set(gcf,'color','w');


%% display the words and columns that are combined in the best bases

% find the levelslengths
[~,LL_rows] = BSfull(GProws,BSrows);
[~,LL_cols] = BSfull(GPcols,BScols);

% print the rows (words) that are combined
IX_rows = find(LL_rows > 1);
fprintf('\n\nCombined words (rows)\n');
fprintf('---------------------\n');
count = 1;
for IX = 1:length(IX_rows)
    fprintf('%4d: %s\n',ind_rows(IX_rows(IX)),words(ind_rows(IX_rows(IX)),:));
    count = count+1;
    if count > LL_rows(IX_rows(IX))
        fprintf('\n');
        count = 1;
    end
end

% print the columns (papers) that are combined
IX_cols = find(LL_cols > 1);
fprintf('\n\nCombined documents (columns)\n');
fprintf('----------------------------\n');
count = 1;
for IX = 1:length(IX_cols)
    fprintf('%4d: %s\n',ind_cols(IX_cols(IX)),articles{ind_cols(IX_cols(IX))});
    count = count+1;
    if count > LL_cols(IX_cols(IX))
        fprintf('\n');
        count = 1;
    end
end