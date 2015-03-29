function [dmatrix,GP] = HGLET_Analysis(G,GP,~,~)
% For a GraphSig object 'G', generate the matrix of HGLET expansion
% coefficients
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

if ~exist('GP','var') || ~isa(GP,'GraphPart') % ~strcmp(class(GP),'GraphPart')
    GP = PartitionTreeFiedler(G);
end

[ind,rs,~,~,~,~,~,method] = ExtractData(GP);
clear GP

N = G.length;
[~,jmax] = size(rs);

% allocate space for the expansion coefficients
[~,fcols] = size(G.f);
dmatrix = zeros(N,jmax,fcols);
dmatrix(:,jmax,:) = G.f(ind,:);


%% 1a. Perform the transform ==> eigenvectors of L

if nargin < 3
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
                dmatrix(rs1,j,:) = G.f(ind(rs1),:);
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
                dmatrix(rs1:rs3-1,j,:) = vec' * G.f(indrs,:);

            end
        end
    end


%% 1b. Perform the transform ==> eigenvectors of L_rw

elseif nargin == 3
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
                dmatrix(rs1,j,:) = G.f(ind(rs1),:);
            elseif n > 1
                indrs = ind(rs1:rs3-1);

                % compute the eigenvectors
                if min(sum(G.W(indrs,indrs),2)) > 10^3*eps
                    useLrw = true;
                    
                    %%% eigenvectors of L_rw ==> svd(L_sym)
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
                
                % obtain the expansion coefficients
                if useLrw
                    dmatrix(rs1:rs3-1,j,:) = vec' * bsxfun(@times, ((full(sum(G.W(indrs,indrs),2))).^0.5), G.f(indrs,:) );
                else
                    dmatrix(rs1:rs3-1,j,:) = vec' * G.f(indrs,:);
                end
                
            end
        end
    end


%% 1c. Perform the transform ==> eigenvectors of L_sym

else
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
                dmatrix(rs1,j,:) = G.f(ind(rs1),:);
            elseif n > 1
                indrs = ind(rs1:rs3-1);

                % compute the eigenvectors
                if min(sum(G.W(indrs,indrs),2)) > 10^3*eps
                    
                    %%% eigenvectors of L_sym ==> svd(L_sym)
                    [vec,~,~] = svd( full( bsxfun(@times, bsxfun(@times, full(sum(G.W(indrs,indrs),2)).^(-0.5), diag(sum(G.W(indrs,indrs))) - G.W(indrs,indrs)), (full(sum(G.W(indrs,indrs),1))).^(-0.5) ) ) );
                    vec = vec(:,end:-1:1);
                    
                else
                    
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

                % obtain the expansion coefficients
                dmatrix(rs1:rs3-1,j,:) = vec' * G.f(indrs,:);

            end
        end
    end
end


if nargout > 1
    GP = GraphPart(ind,rs,[],[],[],[],[],method);
end


end