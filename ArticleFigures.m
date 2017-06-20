function ArticleFigures
% Generate selected figures from the following publications:
% 
% 1. J. Irion and N. Saito, "Hierarchical Graph Laplacian Eigen
% Transforms," Japan SIAM Letters, vol. 6, pp. 21-24, 2014.
%
% 2. J. Irion and N. Saito, "The Generalized Haar-Walsh Transform," Proc.
% 2014 IEEE Workshop on Statistical Signal Processing, pp. 472-475, 2014.
% 
% 3. J. Irion and N. Saito, "Applied and Computational Harmonic Analysis on
% Graphs and Networks," Wavelets and Sparsity XVI, (M. Papadakis, V. K. 
% Goyal, D. Van De Ville, eds.), Proc. SPIE 9597, Paper #95971F, 
% Invited paper, 2015.
%
% 4. J. Irion and N. Saito, "Efficient Approximation and Denoising of Graph
% Signals Using the Multiscale Basis Dictionaries," to appear in IEEE
% Transactions on Signal and Information Processing over Networks, 2017.
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



choice = menu('List of Reproducible Figures', 'HGLET Figure 1', ...
    'GHWT Figure 1', 'GHWT Figure 2', 'GHWT Figure 3', 'GHWT Figure 4', ...
    'ACHA Figure 3', 'ACHA Figure 4', 'ACHA Figure 5', ...
    'TSIPN Figure 2', 'TSIPN Figure 3', 'TSIPN Figure 4', ...
    'TSIPN Figure 5', 'TSIPN Figure 6', 'TSIPN Figure 7', ...
    'TSIPN Figure 8', 'TSIPN Table 1', 'Exit');

switch choice
    case 1
        HGLET_Figure1
        ArticleFigures
    case 2
        GHWT_Figure1
        ArticleFigures
    case 3
        GHWT_Figure2
        ArticleFigures
    case 4
        GHWT_Figure3
        ArticleFigures
    case 5
        GHWT_Figure4
        ArticleFigures
    case 6
        ACHA_Figure3
        ArticleFigures
    case 7
        ACHA_Figure4
        ArticleFigures
    case 8
        ACHA_Figure5
        ArticleFigures
    case 9
        TSIPN_Figure2
        ArticleFigures
    case 10
        TSIPN_Figure3
        ArticleFigures
    case 11
        TSIPN_Figure4
        ArticleFigures
    case 12
        TSIPN_Figure5
        ArticleFigures
    case 13
        TSIPN_Figure6
        ArticleFigures
    case 14
        TSIPN_Figure7
        ArticleFigures
    case 15
        TSIPN_Figure8
        ArticleFigures
    case 16
        TSIPN_Table1
        ArticleFigures
    case 17
        close all
end


end
