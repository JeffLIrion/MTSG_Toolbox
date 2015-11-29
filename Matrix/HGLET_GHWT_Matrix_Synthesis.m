function matrix = HGLET_GHWT_Matrix_Synthesis(dvec,GProws,GPcols,BSrows,BScols,Wrows,Wcols,transrows,transcols)
% Given a vector of HGLET & GHWT expansion coefficients, info about the 
% graph partitioning, and the choice of basis and corresponding transforms,
% reconstruct the signal
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
%   transrows   specifies which transform was used for that portion of the
%               signal (rows): 
%                   00 = HGLET with L
%                   01 = HGLET with Lrw
%                   10 = HGLET with Lsym
%                   11 = GHWT
%   transcols   specifies which transform was used for that portion of the
%               signal (columns)
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

% fill out trans
[~,~,transfullrows] = BSfull(GProws,BSrows,transrows);
[~,~,transfullcols] = BSfull(GPcols,BScols,transcols);


% decompose dvec into GHWT and HGLET components (columns)
dvecH    = bsxfun(@times, dvec', ~transfullcols(:,1) .* ~transfullcols(:,2));
dvecHrw  = bsxfun(@times, dvec', ~transfullcols(:,1) .*  transfullcols(:,2));
dvecHsym = bsxfun(@times, dvec',  transfullcols(:,1) .* ~transfullcols(:,2));
dvecG    = bsxfun(@times, dvec',  transfullcols(:,1) .*  transfullcols(:,2));

% synthesize the row expansion coefficients using the transforms separately
dvecH    = HGLET_Synthesis(dvecH,GPcols,BScols,Gcols);
dvecHrw  = HGLET_Synthesis(dvecHrw,GPcols,BScols,Gcols,1);
dvecHsym = HGLET_Synthesis(dvecHsym,GPcols,BScols,Gcols,1,1);
dvecG    =  GHWT_Synthesis(dvecG,GPcols,BScols);

% combine the row expansion coefficients
dvec = dvecH + dvecHrw + dvecHsym + dvecG;


% decompose dvec into GHWT and HGLET components (rows)
dvecH    = bsxfun(@times, dvec', ~transfullrows(:,1) .* ~transfullrows(:,2));
dvecHrw  = bsxfun(@times, dvec', ~transfullrows(:,1) .*  transfullrows(:,2));
dvecHsym = bsxfun(@times, dvec',  transfullrows(:,1) .* ~transfullrows(:,2));
dvecG    = bsxfun(@times, dvec',  transfullrows(:,1) .*  transfullrows(:,2));

% synthesize the original matrix using the transforms separately
dvecH    = HGLET_Synthesis(dvecH,GProws,BSrows,Grows);
dvecHrw  = HGLET_Synthesis(dvecHrw,GProws,BSrows,Grows,1);
dvecHsym = HGLET_Synthesis(dvecHsym,GProws,BSrows,Grows,1,1);
dvecG    =  GHWT_Synthesis(dvecG,GProws,BSrows);

% combine the row matrix components
matrix = dvecH + dvecHrw + dvecHsym + dvecG;


end