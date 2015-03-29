function Gout = HGLET_BasisVector(G,GP,drow,dcol,~,~)
% Generate the basis vector corresponding to entry dmatrix(drow,dcol).  
%
% Inputs
%   G       a GraphSig object
%   GP      a GraphPart object
%   drow    the row with a "1" entry
%   dcol    the column with a "1" entry
%   ~       if 5 inputs are given, use the eigenvectors of L_rw as the
%           bases
%   ~       if 6 inputs are given, use the eigenvectors of L_sym as the
%           bases
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
if nargin == 4
    [~,Gout] = HGLET_Synthesis(dvec,GP,BS,G);
elseif nargin == 5
    [~,Gout] = HGLET_Synthesis(dvec,GP,BS,G,1);
else
    [~,Gout] = HGLET_Synthesis(dvec,GP,BS,G,1,1);
end

% figure out the j, k, and l indices of the basis vector
[j,k,l] = HGLET_jkl(GP,drow,dcol);

% name the basis vector
if nargin == 4
    Gout = EditName(Gout, sprintf('HGLET (L) Basis Vector (j = %d, k = %d, l = %d)',j,k,l));
elseif nargin == 5
    Gout = EditName(Gout, sprintf('HGLET (Lrw) Basis Vector (j = %d, k = %d, l = %d)',j,k,l));
else
    Gout = EditName(Gout, sprintf('HGLET (Lsym) Basis Vector (j = %d, k = %d, l = %d)',j,k,l));
end

% symmetrize the plot
[~,~,~,~,plotspecs] = ExtractData(Gout);
Gout = EditPlotSpecs(Gout,strcat(plotspecs,'symm'));


end