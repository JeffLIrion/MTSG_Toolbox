function ArticleFigures
% Generate selected figures from the following articles:
% 
% 1. J. Irion and N. Saito, "Hierarchical Graph Laplacian Eigen
% Transforms," Japan SIAM Letters, vol. 6, pp. 21-24, 2014.
%
% 2. J. Irion and N. Saito, "The Generalized Haar-Walsh Transform," Proc. 
% 2014 IEEE Statistical Signal Processing Workshop, pp. 488-491, 2014.
% 
% 3. J. Irion and N. Saito, "Applied and Computational Harmonic Analysis on
% Graphs and Networks," Wavelets and Sparsity XVI, (M. Papadakis, V. K. 
% Goyal, D. Van De Ville, eds.), Proc. SPIE 9597, Paper #95971F, 
% Invited paper, 2015.
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



choice = menu('List of Reproducible Figures', 'HGLET Figure 1', ...
    'GHWT Figure 1', 'GHWT Figure 2', 'GHWT Figure 3', 'GHWT Figure 4', ...
    'ACHA Figure 3', 'ACHA Figure 4', 'ACHA Figure 5', 'Exit');

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
        close all
end


end