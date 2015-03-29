function Gout = GHWT_BasisVector(G,GP,drow,dcol)
% Generate the basis vector corresponding to entry dmatrix(drow,dcol).  
%
% Inputs
%   G       a GraphSig object
%   GP      a GraphPart object
%   drow    the row with a "1" entry
%   dcol    the column with a "1" entry
%
% Outputs
%   Gout    the basis vector corresponding to dmatrix(drow,dcol)
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% get constants
N = G.length;
[~,rs] = ExtractData(GP);
[~,jmax] = size(rs);

% fill in a 1 in the matrix "dmatrix"
dmatrix = zeros(N,jmax);
dmatrix(drow,dcol) = 1;

% convert dmatrix to dvec
[dvec,BS] = dmatrix2dvec(dmatrix,GP);

% generate the basis vector
[~,Gout] = GHWT_Synthesis(dvec,GP,BS,G);

% figure out the j, k, and l indices of the basis vector
[j,k,l] = GHWT_jkl(GP,drow,dcol);

% name the basis vector
Gout = EditName(Gout, sprintf('GHWT Basis Vector (j = %d, k = %d, l = %d)',j,k,l));

% symmetrize the plot
[~,~,~,~,plotspecs] = ExtractData(Gout);
Gout = EditPlotSpecs(Gout,strcat(plotspecs,'symm'));


end