function fig = Matrix_BasisVisual(matrix,GProws,GPcols,BSrows,BScols,rowcolor,colcolor,lw_rows,lw_cols)
% Display an HGLET/GHWT basis for a matrix
%
% Input
%   matrix      the matrix being analyzed
%   GProws      the recursive partitioning on the rows
%   GPcols      the recursive partitioning on the columns
%   BSrows      the basis on the rows (optional)
%   BScols      the basis on the columns (optional)
%   rowcolor    the color to draw the row partitioning lines
%   colcolor    the color to draw the column partitioning lines
%   lw_rows     the LineWidth for the row basis partitions
%   lw_cols     the LineWidth for the column basis partitions

%
% Output
%   fig         the matrix basis visualization
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% If the row and column bases are specified, make sure that both
% are from their respective coarse-to-fine dictionary
if exist('BSrows','var') && isa(BSrows,'BasisSpec') && exist('BScols','var') && isa(BScols,'BasisSpec')
    [~,~,c2f_rows] = ExtractData(BSrows);
    [~,~,c2f_cols] = ExtractData(BScols);
    if ~c2f_rows || ~c2f_cols
        if ~c2f_rows && ~c2f_cols
            fprintf('\n\nThe row and column bases are both from their respective fine-to-coarse dictionary.');
        elseif ~c2f_rows
            fprintf('\n\nThe row basis is from its respective fine-to-coarse dictionary.');
        else
            fprintf('\n\nThe column basis is from its respective fine-to-coarse dictionary.');
        end
        fprintf('\nMatrix_BasisVisual requires the bases to be from the coarse-to-fine dictionary.');
        fprintf('\n\nYou can plot the bases separately using BasisVisual.m.');
        fprintf('\nExiting now.\n\n\n');
        fig = 0;
        return
    end
end

% constants
[rows,cols] = size(matrix);

% extract row & column ordering data
ind_rows = ExtractData(GProws);
ind_cols = ExtractData(GPcols);

% the partitioning trees
[W_GProws,xy_GProws] = PartitionTreeDisplay(GProws);
close(gcf);
[W_GPcols,xy_GPcols] = PartitionTreeDisplay(GPcols);
close(gcf);

% rotate the row recursive partitioning tree
xy_GProws = ([0,1;-1,0]*xy_GProws')';
xy_GProws(:,2) = -xy_GProws(:,2);

% position and stretch the recursive partitioning trees
Nmean = sqrt(rows*cols);
xy_GProws(:,1) =  0.1*Nmean*(xy_GProws(:,1))/log2(Nmean)+cols+0.5;
xy_GPcols(:,2) = -0.1*Nmean*(xy_GPcols(:,2))/log2(Nmean)+0.5;

% plot the matrix
fig = figure;
imagesc(matrix(ind_rows,ind_cols));
xlim([0, max(xy_GProws(:,1))+1]);
ylim([min(xy_GPcols(:,2))-1, rows+1]);
axis equal;
axis off;
hold on;

% plot the partition trees and the row and column bases
if exist('BSrows','var') && isa(BSrows,'BasisSpec') && exist('BScols','var') && isa(BScols,'BasisSpec')
    % plot the partition trees
    gplot(W_GProws,xy_GProws,'-b');
    gplot(W_GPcols,xy_GPcols,'-b');
    
    % specify colors, if necessary
    if ~exist('rowcolor','var')
        rowcolor = 'r';
    end
    if ~exist('colcolor','var')
        colcolor = 'c';
    end
    
    % specify LineWidths, if necessary
    if ~exist('lw_rows','var')
        lw_rows = 0.5;
    end
    if ~exist('lw_cols','var')
        lw_cols = 0.5;
    end
    
    % extract the levlengths data
    [~,levlengthsR] = ExtractData(BSrows); levlengthsR = cumsum(double(levlengthsR(1:end-1)))+0.5;
    [~,levlengthsC] = ExtractData(BScols); levlengthsC = cumsum(double(levlengthsC(1:end-1)))+0.5;

    % plot the row partitions
    x = repmat([0;cols+1],1,length(levlengthsR));
    plot(x,[levlengthsR';levlengthsR'],'-','Color',rowcolor,'LineWidth',lw_rows);

    % plot the column partitions
    y = repmat([0;rows+1],1,length(levlengthsC));
    plot([levlengthsC';levlengthsC'],y,'-','Color',colcolor,'LineWidth',lw_cols);
    
    % use a gray colormap
    colormap(gray);
    
else
    % plot the partition trees
    gplot(W_GProws,xy_GProws,'-k');
    gplot(W_GPcols,xy_GPcols,'-k');
    
    % use the jet colormap if BSrows and BScols are not provided
    colormap(jet);
end

set(gcf,'Color','w');


end