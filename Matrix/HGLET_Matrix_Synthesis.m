function matrix = HGLET_Matrix_Synthesis(dvec,GProws,GPcols,BSrows,BScols,Wrows,Wcols,~,~)
% Given a vector of HGLET expansion coefficients and info about the graph 
% partitioning and the choice of basis, reconstruct the signal
%
% Input
%   dvec        the expansion coefficients corresponding to the chosen
%               basis
%   GProws      the recursive partitioning on the rows
%   GPcols      the recursive partitioning on the columns
%   BSrows      the best basis on the rows
%   BScols      the best basis on the columns
%   Wrows       the weight matrix on the rows
%   Wcols       the weight matrix on the columns
%   ~           if 8 inputs are given, use the eigenvectors of L_rw as the
%               bases  ==>  U_rw = D^(0.5) * U_sym
%   ~           if 9 inputs are given, use the eigenvectors of L_sym as the
%               bases  ==>  U_sym = D^(-0.5) * U_rw
%
% Output
%   matrix      the reconstructed matrix
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% generate GraphSig objects for the matrix
Grows = GraphSig(Wrows,[],[]);
Gcols = GraphSig(Wcols,[],[]);

% reconstruct the columns
if nargin < 8
    dvec = HGLET_Synthesis(dvec',GPcols,BScols,Gcols);
elseif nargin == 8
    dvec = HGLET_Synthesis(dvec',GPcols,BScols,Gcols,1);
else
    dvec = HGLET_Synthesis(dvec',GPcols,BScols,Gcols,1,1);
end

% reconstruct the rows
if nargin < 8
    matrix = HGLET_Synthesis(dvec',GProws,BSrows,Grows);
elseif nargin == 8
    matrix = HGLET_Synthesis(dvec',GProws,BSrows,Grows,1);
else
    matrix = HGLET_Synthesis(dvec',GProws,BSrows,Grows,1,1);
end


end