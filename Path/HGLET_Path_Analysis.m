function [dmatrix,GP] = HGLET_Path_Analysis(G,GP,~,~)
% For a GraphSig object 'G' corresponding to an unweighted path graph, 
% generate the matrix of HGLET expansion coefficients
%
% DCT Info:
%  http://vadkudr.org/Algorithms/DTT/DCTI_DCTII/DCTI_DCTII.html#5
%  "Adapted Wavelet Analysis: From Theory to Software" by M.V. Wickerhauser
%  (page 88)
%
% Input
%   G               a GraphSig object
%   GP              a GraphPart object
%   ~               if a 3rd input is given, use the eigevectors of L_rw as
%                   the bases  ==>  U_rw = D^(0.5) * U_sym
%   ~               if a 4th input is given, use the eigenvectors of L_sym
%                   as the bases  ==>  U_sym = D^(-0.5) * U_rw
%
% Output
%   dmatrix         the matrix of expansion coefficients
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
dmatrix = zeros(N,jmax,fcols);
dmatrix(:,jmax,:) = f(ind,:);


%% 1a. Perform the transform ==> eigenvectors of L

if nargin < 3
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
                dmatrix(rs1,j,:) = f(ind(rs1),:);
            elseif n > 1
                indrs = ind(rs1:rs3-1);

                % obtain the expansion coefficients ==> DCT-II
                dmatrix(rs1:rs3-1,j,:) = dct(f(indrs,:));
            end
        end
    end


%% 1b. Perform the transform ==> eigenvectors of L_rw

elseif nargin == 3
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
                dmatrix(rs1,j,:) = f(ind(rs1),:);
            elseif n > 1
                indrs = ind(rs1:rs3-1);
                
                % obtain the expansion coefficients ==> DCT-I
                D = 2*ones(n,1);
                D([1,n]) = 1;
                dmatrix(rs1:rs3-1,j,:) = dct1(bsxfun(@times, D.^0.5, f(indrs,:)));
            end
        end
    end


%% 1c. Perform the transform ==> eigenvectors of L_sym

else
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
                dmatrix(rs1,j,:) = f(ind(rs1),:);
            elseif n > 1
                indrs = ind(rs1:rs3-1);

                % obtain the expansion coefficients ==> DCT-I
                dmatrix(rs1:rs3-1,j,:) = dct1(f(indrs,:));
            end
        end
    end
end


if nargout > 1
    GP = GraphPart(ind,rs,[],[],[],[],[],method);
end


end
