function [f,GS] = HGLET_Synthesis(dvec,GP,BS,G,~,~)
% Given a vector of HGLET expansion coefficients and info about the graph 
% partitioning and the choice of basis, reconstruct the signal
%
% Input
%   dvec        the expansion coefficients corresponding to the chosen
%               basis
%   GP          a GraphPart object
%   BS          a BasisSpec object
%   G           a GraphSig object
%   ~           if a 5th input is given, use the eigenvectors of L_rw as
%               the bases  ==>  U_rw = D^(0.5) * U_sym
%   ~           if a 6th input is given, use the eigenvectors of L_sym as
%               the bases  ==>  U_sym = D^(-0.5) * U_rw
%
% Output
%   f           the reconstructed signal
%   GS          the reconstructed GraphSig object
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminaries

% constants
[~,jmax] = size(GP.rs);

% fill in the appropriate entries of dmatrix
dmatrix = dvec2dmatrix(dvec,GP,BS);
f = squeeze(dmatrix(:,jmax,:));

W = ExtractData(G);


%% 1a. Perform the synthesis transform ==> eigenvectors of L

if nargin < 5
    for j = jmax:-1:1
        regioncount = nnz(GP.rs(:,j))-1;
        for r = 1:regioncount
            % the index that marks the start of the region
            rs1 = GP.rs(r,j);

            % the index that is one after the end of the region
            rs3 = GP.rs(r+1,j);

            % the number of points in the current region
            n = double(rs3 - rs1);

            % only proceed forward if coefficients do not exist
            if (j == jmax || nnz(dmatrix(rs1:rs3-1,j+1,:)) == 0) && nnz(dmatrix(rs1:rs3-1,j,:)) > 0

                if n == 1
                    f(rs1,:) = dmatrix(rs1,j,:);
                elseif n > 1
                    indrs = GP.ind(rs1:rs3-1);

                    % compute the eigenvectors of L ==> svd(L)
                    [vec,~,~] = svd(full( diag(sum(W(indrs,indrs))) - W(indrs,indrs) ));
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

                    % reconstruct the signal
                    f(rs1:rs3-1,:) = vec * squeeze(dmatrix(rs1:rs3-1,j,:));

                end

            end
        end
    end


%% 1b. Perform the synthesis transform ==> eigenvectors of L_rw

elseif nargin == 5
    for j = jmax:-1:1
        regioncount = nnz(GP.rs(:,j))-1;
        for r = 1:regioncount
            % the index that marks the start of the region
            rs1 = GP.rs(r,j);

            % the index that is one after the end of the region
            rs3 = GP.rs(r+1,j);

            % the number of points in the current region
            n = double(rs3 - rs1);

            % only proceed forward if coefficients do not exist
            if (j == jmax || nnz(dmatrix(rs1:rs3-1,j+1,:))) == 0 && nnz(dmatrix(rs1:rs3-1,j,:)) > 0

                if n == 1
                    f(rs1,:) = dmatrix(rs1,j,:);
                elseif n > 1
                    indrs = GP.ind(rs1:rs3-1);

                    % compute the eigenvectors
                    if min(sum(W(indrs,indrs),2)) > 10^3*eps
                        useLrw = true;
                        
                        %%% eigenvectors of L_rw ==> svd(L_sym)
                        [vec,~,~] = svd( full( bsxfun(@times, bsxfun(@times, full(sum(W(indrs,indrs),2)).^(-0.5), diag(sum(W(indrs,indrs))) - W(indrs,indrs)), (full(sum(W(indrs,indrs),1))).^(-0.5) ) ) );
                        vec = vec(:,end:-1:1);
                        
                    else
                        useLrw = false;
                        
                        %%% eigenvectors of L ==> svd(L)
                        [vec,~,~] = svd(full( diag(sum(W(indrs,indrs))) - W(indrs,indrs) ));
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
                    
                    % reconstruct the signal
                    if useLrw
                        f(rs1:rs3-1,:) = bsxfun(@times, (full(sum(W(indrs,indrs),2))).^(-0.5), vec * squeeze(dmatrix(rs1:rs3-1,j,:)) );
                    else
                        f(rs1:rs3-1,:) = vec * squeeze(dmatrix(rs1:rs3-1,j,:));
                    end

                end

            end
        end
    end


%% 1c. Perform the synthesis transform ==> eigenvectors of L_sym

else
    for j = jmax:-1:1
        regioncount = nnz(GP.rs(:,j))-1;
        for r = 1:regioncount
            % the index that marks the start of the region
            rs1 = GP.rs(r,j);

            % the index that is one after the end of the region
            rs3 = GP.rs(r+1,j);

            % the number of points in the current region
            n = double(rs3 - rs1);

            % only proceed forward if expansion coefficients exist and
            % signal values do not
            if (j == jmax || nnz(dmatrix(rs1:rs3-1,j+1,:)) == 0) && nnz(dmatrix(rs1:rs3-1,j,:)) > 0

                if n == 1
                    f(rs1,:) = dmatrix(rs1,j,:);
                elseif n > 1
                    indrs = GP.ind(rs1:rs3-1);

                    % compute the eigenvectors
                    if min(sum(W(indrs,indrs),2)) > 10^3*eps
                        
                        %%% eigenvectors of L_sym ==> svd(L_sym)
                        [vec,~,~] = svd( full( bsxfun(@times, bsxfun(@times, full(sum(W(indrs,indrs),2)).^(-0.5), diag(sum(W(indrs,indrs))) - W(indrs,indrs)), (full(sum(W(indrs,indrs),1))).^(-0.5) ) ) );
                        vec = vec(:,end:-1:1);
                        
                    else
                        
                        %%% eigenvectors of L ==> svd(L)
                        [vec,~,~] = svd(full( diag(sum(W(indrs,indrs))) - W(indrs,indrs) ));
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

                    % reconstruct the signal
                    f(rs1:rs3-1,:) = vec * squeeze(dmatrix(rs1:rs3-1,j,:));

                end
            end
        end
    end
end



% put the reconstructed values in the correct order
f(GP.ind,:) = f;


% create a GraphSig object with the reconstructed data
if nargout == 2
    if exist('G','var') == 1 && isa(G,'GraphSig')
        GS = ReplaceData(G,f);
    else
        GS = 0;
    end
end


end