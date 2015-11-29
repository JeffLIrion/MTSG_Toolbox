function [dmatrixH,dmatrixHrw,dmatrixHsym,GP] = HGLET_Path_Analysis_All(G,GP)
% For a GraphSig object 'G' corresponding to an unweighted path graph, 
% generate the 3 matrices of HGLET expansion coefficients corresponding to
% the eigenvectors of L, Lrw, and Lsym
%
% DCT Info:
%  http://vadkudr.org/Algorithms/DTT/DCTI_DCTII/DCTI_DCTII.html#5
%  "Adapted Wavelet Analysis: From Theory to Software" by M.V. Wickerhauser
%  (page 88)
%
% Input
%   G               a GraphSig object
%   GP              a GraphPart object
%
% Output
%   dmatrixH        the matrix of expansion coefficients for L
%   dmatrixHrw      the matrix of expansion coefficients for Lrw
%   dmatrixHsym     the matrix of expansion coefficients for Lsym
%   GP              a GraphPart object
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminaries

if ~exist('GP','var') || ~isa(GP,'GraphPart')
    GP = PartitionTreeFiedler(G);
end

[~,~,f] = ExtractData(G);
[ind,rs,~,~,~,~,~,method] = ExtractData(GP);
clear GP

N = length(ind);
[~,jmax] = size(rs);

% allocate space for the expansion coefficients
[~,fcols] = size(f);
dmatrixH = zeros(N,jmax,fcols);
dmatrixH(:,jmax,:) = f(ind,:);
dmatrixHrw = dmatrixH;
dmatrixHsym = dmatrixH;


%% 1a. Perform the transform ==> eigenvectors of L

for j = jmax:-1:1
    regioncount = nnz(rs(:,j))-1;
    for r = 1:regioncount
        % the index that marks the start of the region
        rs1 = rs(r,j);

        % the index that is one after the end of the region
        rs3 = rs(r+1,j);

        % the number of points in the current region
        n = double(rs3 - rs1);

        if n == 1
            dmatrixH(rs1,j,:) = f(ind(rs1),:);
        elseif n > 1
            indrs = ind(rs1:rs3-1);

            % obtain the expansion coefficients ==> DCT-II
            dmatrixH(rs1:rs3-1,j,:) = dct(f(indrs,:));
        end
    end
end


%% 1b. Perform the transform ==> eigenvectors of L_rw and Lsym

for j = jmax:-1:1
    regioncount = nnz(rs(:,j))-1;
    for r = 1:regioncount
        % the index that marks the start of the region
        rs1 = rs(r,j);

        % the index that is one after the end of the region
        rs3 = rs(r+1,j);

        % the number of points in the current region
        n = double(rs3 - rs1);

        if n == 1
            dmatrixHrw(rs1,j,:) = f(ind(rs1),:);
            dmatrixHsym(rs1,j,:) = f(ind(rs1),:);
        elseif n > 1
            indrs = ind(rs1:rs3-1);
            
            % obtain the expansion coefficients ==> DCT-I
            D = 2*ones(n,1);
            D([1,n]) = 1;
            dmatrixHrw(rs1:rs3-1,j,:) = dct1(bsxfun(@times, D.^0.5, f(indrs,:)));
            dmatrixHsym(rs1:rs3-1,j,:) = dct1(f(indrs,:));
        end
    end
end



if nargout > 1
    GP = GraphPart(ind,rs,[],[],[],[],[],method);
end


end