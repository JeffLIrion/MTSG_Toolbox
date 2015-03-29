function [Vis,GVis] = GHWT_BasisVisual(G,GP,BS,dvec)
% Display a GHWT basis
%
% Input
%   G           a GraphSig object
%   GP          a GraphPart object
%   BS          a BasisSpec object
%   dvec        the expansion coefficients corresponding to the chosen
%               basis
%
% Output
%   Vis         a visualization of the specified basis ==> imagesc(Vis)
%   GVis        a GraphSig object that displays the specified basis
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% extract data and figure out which dictionaries I'm dealing with
[ind,rs,tag,~,~,tagf2c] = ExtractData(GP);
[N,jmax] = size(rs);
N = N-1;

if isempty(tagf2c)
    GP = GHWT_Info(GP);
    [~,~,tag,~,~,tagf2c] = ExtractData(GP);
end
[~,~,BSc2f] = ExtractData(BS);

% extract levlist and levlengths info
[levlist,levlengths] = ExtractData(BS);
if isempty(levlengths)
    BS = levlist2levlengths(GP,BS);
    [levlist,levlengths] = ExtractData(BS);
end

% allocate spece
Vis = zeros(jmax,G.length);


% coarse-to-fine
if BSc2f
    % fill in the matrix
    indr0 = 1;
    for row = 1:length(levlist)
        indr = indr0:indr0+levlengths(row)-1;
        if exist('dvec','var') == 1
            Vis(levlist(row),indr) = abs(dvec(indr));
        else
            Vis(levlist(row),indr) = 1 + log2(double(tag(indr,levlist(row))+1));
        end
        indr0 = indr0+levlengths(row);
    end
    
    % generate the GraphSig object
    if nargout == 2
        W = sparse(G.length,G.length);

        % modify the weight matrix
        indr0 = 1;
        for row = 1:length(levlist)
            indr = indr0:indr0+levlengths(row)-1;
            inds = ind(indr);
            W(inds,inds) = G.W(inds,inds);
            indr0 = indr0+levlengths(row);
        end

        G.W = W;
        if exist('dvec','var') == 1
            G.f = zeros(G.length,1);
            G.f(ind) = sum(abs(dvec),2);
        else
            G.f = zeros(G.length,1);
            G.f(ind) = BSfull(GP,BS);
        end

        GVis = G;
    end
    
    
% fine-to-coarse
else
    
    % fill in the matrix
    indr0 = 1;
    for row = 1:length(levlist)
        indr = indr0:indr0+levlengths(row)-1;
        if exist('dvec','var') == 1
            Vis(levlist(row),indr) = abs(dvec(indr));
        else
            Vis(levlist(row),indr) = ceil(log2(double(tagf2c(indr0+levlengths(row)-1,levlist(row))+1)));
        end
        indr0 = indr0+levlengths(row);
    end
    
    % generate the GraphSig object
    if nargout == 2
        if ~exist('dvec','var')
            GVis = G;
        else
            % see which coefficient contributes the most to each node
            
            % fill in the appropriate entries of dmatrix
            dmatrix = dvec2dmatrix(dvec,GP,BS);
            
            % constants
            [~,fcols] = size(dvec);
            
            % extract info
            [ind,tag,~,~,rsf2c,tagf2c,compinfof2c] = ExtractData(GP);
            
            % obtain the fine-to-coarse indexing scheme
            [~,~,IX] = FineToCoarse(GP,zeros(N,jmax));
            
            % perform the synthesis transform
            for j = jmax-1:-1:1
                regioncount = nnz(rsf2c(:,j))-1;
                for r = 1:regioncount
                    % the index that marks the start of the first subregion
                    rs1 = rsf2c(r,j);

                    % the index that is one after the end of the second subregion
                    rs3 = rsf2c(r+1,j);

                    % only proceed forward if coefficients do not exist
                    if nnz(dmatrix(rs1:rs3-1,j,:)) == 0 && nnz(dmatrix(rs1:rs3-1,j+1,:)) > 0

                        % one subregion ==> copy the coefficients
                        if tagf2c(rs1,j+1) == tagf2c(rs3-1,j+1)
                            dmatrix(rs1:rs3-1,j,:) = dmatrix(rs1:rs3-1,j+1,:);

                        % two subregions ==> compute the coefficients
                        else
                            % the index that marks the start of the second subregion
                            rs2 = rs1+1;
                            while tagf2c(rs1,j+1) == tagf2c(rs2,j+1) && rs2 < N+1
                                rs2 = rs2+1;
                            end

                            % the current coefficient in the parent region with which we  are dealing
                            parent = rs1;

                            % the current coefficients in the child regions with which we are dealing
                            child1 = rs1;
                            child2 = rs2;


                            %%% SCALING COEFFICIENTS
                            if tagf2c(rs1,j) == 0
                                while parent < rs3
                                    % 1 coefficient formed from 1 coefficient
                                    if compinfof2c(child1,j+1) == 0
                                        dmatrix(parent,j,:) = dmatrix(child1,j+1,:);
                                        parent = parent+1;
                                        child1 = child1+1;

                                    % 2 coefficients formed from 2 coefficients
                                    else
                                        % the number of points in the first subregion
                                        n1 = double(compinfof2c(child1,j+1));

                                        % the number of points in the second subregion
                                        n2 = double(compinfof2c(child2,j+1));

                                        % the number of points in the parent region
                                        n = n1+n2;

                                        dmatrix(parent,  j,:) = ( sqrt(n1)*dmatrix(child1,j+1,:) + sqrt(n2)*dmatrix(child2,j+1,:) )/sqrt(n);
                                        dmatrix(parent+1,j,:) = ( sqrt(n2)*dmatrix(child1,j+1,:) - sqrt(n1)*dmatrix(child2,j+1,:) )/sqrt(n);
                                        parent = parent+2;
                                        child1 = child1+1;
                                        child2 = child2+1;
                                    end
                                end

                            %%% HAAR-LIKE & WALSH-LIKE COEFFICIENTS
                            else
                                while parent < rs3
                                    % 1 coefficient formed from 1 coefficient
                                    if compinfof2c(child1,j+1) == 0
                                        dmatrix(parent,j,:) = dmatrix(child1,j+1,:);
                                        parent = parent+1;
                                        child1 = child1+1;

                                    % 2 coefficients formed from 2 coefficients
                                    else
                                        dmatrix(parent  ,j,:) = ( dmatrix(child1,j+1,:) + dmatrix(child2,j+1,:) )/sqrt(2);
                                        dmatrix(parent+1,j,:) = ( dmatrix(child1,j+1,:) - dmatrix(child2,j+1,:) )/sqrt(2);
                                        parent = parent+2;
                                        child1 = child1+1;
                                        child2 = child2+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % put dmatrix in coarse-to-fine order
            for j = 1:jmax
                dmatrix(IX(:,j),j,:) = dmatrix(:,j,:);
            end
                
            % find the basis vectors that contribute the most to each node
            f = zeros(N,fcols);
            
            for col = 1:fcols
                for row = 1:N
                    % level
%                     [~,f(row,col)] = max(abs(dmatrix(row,:,col)));
                    
                    % log(tag)
                    [~,IX] = max(abs(dmatrix(row,:,col)));
                    f(row,col) = log2(double(tag(row,IX)+1));
                end
            end
            
            f = f-1;
            
            % put the reconstructed values in the correct order
            f(ind,:) = f;
            
            GVis = ReplaceData(G,f);
    
        end
    end
end


end