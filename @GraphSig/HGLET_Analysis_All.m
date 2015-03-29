function [dmatrixH,dmatrixHrw,dmatrixHsym,GP] = HGLET_Analysis_All(G,GP)
% For a GraphSig object 'G', generate the 3 matrices of HGLET expansion 
% coefficients corresponding to the eigenvectors of L, Lrw, and Lsym
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

[ind,rs,~,~,~,~,~,method] = ExtractData(GP);
clear GP

N = G.length;
[~,jmax] = size(rs);

% allocate space for the expansion coefficients
[~,fcols] = size(G.f);
dmatrixH = zeros(N,jmax,fcols);
dmatrixH(:,jmax,:) = G.f(ind,:);
dmatrixHrw = dmatrixH;
dmatrixHsym = dmatrixH;


%% 1a. Perform the transform ==> eigenvectors of L

for j = jmax-1:-1:1
    regioncount = nnz(rs(:,j))-1;
    for r = 1:regioncount
        % the index that marks the start of the region
        rs1 = rs(r,j);

        % the index that is one after the end of the region
        rs3 = rs(r+1,j);

        % the number of points in the current region
        n = double(rs3 - rs1);

        if n == 1
            dmatrixH(rs1,j,:) = G.f(ind(rs1),:);
        elseif n > 1
            indrs = ind(rs1:rs3-1);

            % compute the eigenvectors of L ==> svd(L)
            [vec,~,~] = svd(full( diag(sum(G.W(indrs,indrs))) - G.W(indrs,indrs) ));
            vec = vec(:,end:-1:1);

            % standardize the eigenvector signs
            for col = 1:n
                row = 1;
                standardized = false;
                while ~standardized
                    if vec(row,col) > 10^3*eps
                        standardized = true;
                    elseif vec(row,col) < -10^3*eps
                        vec(:,col) = -vec(:,col);
                        standardized = true;
                    else
                        row = row+1;
                    end
                end
            end

            % obtain the expansion coefficients
            dmatrixH(rs1:rs3-1,j,:) = vec' * G.f(indrs,:);

        end
    end
end


%% 1b. Perform the transform ==> eigenvectors of L_rw and Lsym

for j = jmax-1:-1:1
    regioncount = nnz(rs(:,j))-1;
    for r = 1:regioncount
        % the index that marks the start of the region
        rs1 = rs(r,j);

        % the index that is one after the end of the region
        rs3 = rs(r+1,j);

        % the number of points in the current region
        n = double(rs3 - rs1);

        if n == 1
            dmatrixHrw(rs1,j,:) = G.f(ind(rs1),:);
            dmatrixHsym(rs1,j,:) = G.f(ind(rs1),:);
        elseif n > 1
            indrs = ind(rs1:rs3-1);

            % compute the eigenvectors
            if min(sum(G.W(indrs,indrs),2)) > 10^3*eps
                useLrw = true;

                %%% eigenvectors of L_rw and L_sym ==> svd(L_sym)
                [vec,~,~] = svd( full( bsxfun(@times, bsxfun(@times, full(sum(G.W(indrs,indrs),2)).^(-0.5), diag(sum(G.W(indrs,indrs))) - G.W(indrs,indrs)), (full(sum(G.W(indrs,indrs),1))).^(-0.5) ) ) );
                vec = vec(:,end:-1:1);

            else
                useLrw = false;

                %%% eigenvectors of L ==> svd(L)
                [vec,~,~] = svd(full( diag(sum(G.W(indrs,indrs))) - G.W(indrs,indrs) ));
                vec = vec(:,end:-1:1);
            end

            % standardize the eigenvector signs
            for col = 1:n
                row = 1;
                standardized = false;
                while ~standardized
                    if vec(row,col) > 10^3*eps
                        standardized = true;
                    elseif vec(row,col) < -10^3*eps
                        vec(:,col) = -vec(:,col);
                        standardized = true;
                    else
                        row = row+1;
                    end
                end
            end

            % obtain the expansion coefficients for L_sym
            dmatrixHsym(rs1:rs3-1,j,:) = vec' * G.f(indrs,:);
                
            % obtain the expansion coefficients for L_rw
            if useLrw
                dmatrixHrw(rs1:rs3-1,j,:) = vec' * bsxfun(@times, ((full(sum(G.W(indrs,indrs),2))).^0.5), G.f(indrs,:) );
            else
                dmatrixHrw(rs1:rs3-1,j,:) = vec' * G.f(indrs,:);
            end

        end
    end
end



if nargout > 1
    GP = GraphPart(ind,rs,[],[],[],[],[],method);
end


end