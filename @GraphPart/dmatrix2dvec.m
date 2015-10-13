function [dvec,BS] = dmatrix2dvec(dmatrix,GP,BS)
% Given a matrix of expansion coefficients, convert it to a vector.
%
% Inputs
%   dmatrix     a matrix of expansion coefficients
%   GP          a GraphPart object
%   BS          a BasisSpec object
%
% Outputs
%   dvec        a vector of expansion coefficients
%   BS          a BasisSpec object
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% constants
[N,jmax,fcols] = size(dmatrix);


% THE BASIS IS SPECIFIED ==> select the corresponding coefficients
if exist('BS','var')
    %% 0. Preliminaries

    % allocate space
    dvec = zeros(N,fcols);

    [levlist,levlengths,BSc2f] = ExtractData(BS,GP);
    
    % put dmatrix in the fine-to-coarse arrangement, if necessary
    if ~BSc2f
        [~,dmatrix] = FineToCoarse(GP,dmatrix);
    end


    %% 1. Make a vector of the matrix entries specified by the BasisSpec object
    n = 1;
    for row = 1:length(levlist)
        dvec(n:n+levlengths(row)-1,:) = dmatrix(n:n+levlengths(row)-1,levlist(row),:);
        n = n+levlengths(row);
    end
    
    
% THE BASIS IS NOT SPECIFIED ==> retain nonzero coefficients and specify the basis
else
    %% 0. Preliminaries
    
    % allocate/initialize
    dvec = dmatrix(:,jmax,:);
    levlist = jmax*ones(N,1,'uint8');
    
    %% 1. Make a vector of the nonzero basis entries
    for j = jmax-1:-1:1
        regioncount = nnz(GP.rs(:,j))-1;
        for r = 1:regioncount
            indr = GP.rs(r,j):GP.rs(r+1,j)-1;
            if nnz(dmatrix(indr,j,:)) > 0
                dvec(indr,:) = dmatrix(indr,j,:);
                levlist(GP.rs(r,j)) = j;
                levlist(GP.rs(r,j)+1:GP.rs(r+1,j)-1) = 0;
            end
        end
    end
    
    % specify the corresponding basis
    levlist( levlist==0 ) = [];

    BS = BasisSpec(levlist,[],true);
    BS = levlist2levlengths(GP,BS);
    
end


end