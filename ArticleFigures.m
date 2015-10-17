function ArticleFigures
% Generate selected figures from "Hierarchical Graph Laplacian Eigen
% Transforms" and "The Generalized Haar-Walsh Transform"
% 
% 1. J. Irion and N. Saito, "Hierarchical graph Laplacian eigen
% transforms," Japan SIAM Letters, vol. 6, pp. 21-24, 2014.
% 2. J. Irion and N. Saito, "The generalized Haar-Walsh transform," Proc. 
% 2014 IEEE Statistical Signal Processing Workshop, pp. 488-491, 2014.
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



choice = menu('List of Reproducible Figures', 'HGLET Figure 1', ...
    'GHWT Figure 1', 'GHWT Figure 2', 'GHWT Figure 3', 'GHWT Figure 4', ...
    'Exit');

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
        close all
end


end