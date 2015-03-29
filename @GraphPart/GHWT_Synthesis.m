function [f,GS] = GHWT_Synthesis(dvec,GP,BS,G)
% Given a vector of GHWT expansion coefficients and info about the graph 
% partitioning and the choice of basis, reconstruct the signal
%
% Input
%   dvec        the expansion coefficients corresponding to the chosen
%               basis
%   GP          a GraphPart object
%   BS          a BasisSpec object
%   G           a GraphSig object (optional)
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



% if necessary, fill in the tag info in GP
if isempty(GP.tag)
    GP = GHWT_Core(GP);
end

% figure out which dictionary is used: coarse-to-fine or fine-to-coarse
[~,~,c2f] = ExtractData(BS);
if ~c2f && isempty(GP.tagf2c)
    GP = FineToCoarse(GP);
end

% constants
N = length(GP.ind);
[~,jmax] = size(GP.rs);
[~,fcols] = size(dvec);

% fill in the appropriate entries of dmatrix
dmatrix = dvec2dmatrix(dvec,GP,BS);


%% Synthesis from the coarse-to-fine dictionary

if c2f
    for j = 1:jmax-1
        regioncount = nnz(GP.rs(:,j))-1;
        for r = 1:regioncount
            % the index that marks the start of the first subregion
            rs1 = GP.rs(r,j);

            % the index that is one after the end of the second subregion
            rs3 = GP.rs(r+1,j);

            % the number of points in the current region
            n = double(rs3 - rs1);

            % only proceed forward if coefficients do not exist
            if nnz(dmatrix(rs1:rs3-1,j+1,:)) == 0 && nnz(dmatrix(rs1:rs3-1,j,:)) > 0

                if n == 1
                    %%% SCALING COEFFICIENT (n == 1)
                    dmatrix(rs1,j+1,:) = dmatrix(rs1,j,:);
                elseif n > 1
                    % the index that marks the start of the second subregion
                    rs2 = rs1+1;
                    while GP.tag(rs2,j+1) ~= 0 && rs2 < N+1
                        rs2 = rs2+1;
                    end

                    % the number of points in the first subregion
                    n1 = double(rs2 - rs1);

                    % the number of points in the second subregion
                    n2 = double(rs3 - rs2);

                    %%% SCALING COEFFICIENTS
                    dmatrix(rs1,j+1,:) = ( sqrt(n1)*dmatrix(rs1,j,:) + sqrt(n2)*dmatrix(rs1+1,j,:) )/sqrt(n);
                    dmatrix(rs2,j+1,:) = ( sqrt(n2)*dmatrix(rs1,j,:) - sqrt(n1)*dmatrix(rs1+1,j,:) )/sqrt(n);

                    %%% HAAR-LIKE & WALSH-LIKE COEFFICIENTS

                    % search through the remaining coefficients in each subregion
                    parent = rs1+2;
                    child1 = rs1+1;
                    child2 = rs2+1;
                    while child1 < rs2 || child2 < rs3
                        % subregion 1 has the smaller tag
                        if child2 == rs3 || (GP.tag(child1,j+1) < GP.tag(child2,j+1) && child1 < rs2)
                            dmatrix(child1,j+1,:) = dmatrix(parent,j,:);
                            child1 = child1+1;
                            parent = parent+1;

                        % subregion 2 has the smaller tag
                        elseif child1 == rs2 || (GP.tag(child2,j+1) < GP.tag(child1,j+1) && child2 < rs3)
                            dmatrix(child2,j+1,:) = dmatrix(parent,j,:);
                            child2 = child2+1;
                            parent = parent+1;

                        % both subregions have the same tag
                        else
                            dmatrix(child1,j+1,:) = ( dmatrix(parent,j,:) + dmatrix(parent+1,j,:) )/sqrt(2);
                            dmatrix(child2,j+1,:) = ( dmatrix(parent,j,:) - dmatrix(parent+1,j,:) )/sqrt(2);
                            child1 = child1+1;
                            child2 = child2+1;
                            parent = parent+2;
                        end
                    end
                end

            end
        end
    end
    if fcols == 1
        f = dmatrix(:,end);
    else
        f = reshape(dmatrix(:,end,:),N,fcols);
    end
    
    
%% Synthesis from the fine-to-coarse dictionary
else
    for j = jmax-1:-1:1
        regioncount = nnz(GP.rsf2c(:,j))-1;
        for r = 1:regioncount
            % the index that marks the start of the first subregion
            rs1 = GP.rsf2c(r,j);

            % the index that is one after the end of the second subregion
            rs3 = GP.rsf2c(r+1,j);

            % only proceed forward if coefficients do not exist
            if nnz(dmatrix(rs1:rs3-1,j,:)) == 0 && nnz(dmatrix(rs1:rs3-1,j+1,:)) > 0

                % one subregion ==> copy the coefficients
                if GP.tagf2c(rs1,j+1) == GP.tagf2c(rs3-1,j+1)
                    dmatrix(rs1:rs3-1,j,:) = dmatrix(rs1:rs3-1,j+1,:);
                    
                % two subregions ==> compute the coefficients
                else
                    % the index that marks the start of the second subregion
                    rs2 = rs1+1;
                    while GP.tagf2c(rs1,j+1) == GP.tagf2c(rs2,j+1) && rs2 < N+1
                        rs2 = rs2+1;
                    end

                    % the current coefficient in the parent region with which we  are dealing
                    parent = rs1;
                    
                    % the current coefficients in the child regions with which we are dealing
                    child1 = rs1;
                    child2 = rs2;
                    
                    
                    %%% SCALING COEFFICIENTS
                    if GP.tagf2c(rs1,j) == 0
                        while parent < rs3
                            % 1 coefficient formed from 1 coefficient
                            if GP.compinfof2c(child1,j+1) == 0
                                dmatrix(parent,j,:) = dmatrix(child1,j+1,:);
                                parent = parent+1;
                                child1 = child1+1;
                                
                            % 2 coefficients formed from 2 coefficients
                            else
                                % the number of points in the first subregion
                                n1 = double(GP.compinfof2c(child1,j+1));

                                % the number of points in the second subregion
                                n2 = double(GP.compinfof2c(child2,j+1));
                                
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
                            if GP.compinfof2c(child1,j+1) == 0
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
    if fcols == 1
        f = dmatrix(:,1);
    else
        f = reshape(dmatrix(:,1,:),N,fcols);
    end
end



% put the reconstructed values in the correct order
f(GP.ind,:) = f;


% create a GraphSig object with the reconstructed data
if nargout == 2
    if exist('G','var') == 1
        GS = ReplaceData(G,f);
    else
        GS = 0;
    end
end


end