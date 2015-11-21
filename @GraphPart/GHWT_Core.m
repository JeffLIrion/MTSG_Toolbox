function [GP,dmatrix] = GHWT_Core(GP,dmatrix)
% Perform the forward GHWT transform, which generates expansion
% coefficients, tag, and compinfo.  
%
% Input
%   GP          a GraphPart object (with or without tag and compinfo data)
%   dmatrix     a matrix of expansion coefficients with only the last
%               column filled in
%
% Output
%   GP          a GraphPart object with tag and compinfo data
%   dmatrix     a matrix of expansion coefficients
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminaries

[N,jmax] = size(GP.rs);
N = N-1;

% tag -- when expressed in binary, tag indicates the sequence of 
% average (0's) and difference (1's) operations used to obtain the given 
% basis vector
sslac = tag_class(jmax);
tag = zeros(N,jmax,sslac);

% allocate space for compinfo (used for synthesis from fine-to-coarse
% dictionary)
%       if tag == 0 && n == 1:                              compinfo = 0
%       if tag == 0 && n >= 2:                              compinfo = n1
%       if tag == 1:                                        compinfo = n2
%       if tag >= 2 && coeff is formed from 2 coeffs:       compinfo = 1
%       if tag >= 2 && coeff is formed from 1 coeff:        compinfo = 0
compinfo = zeros(N,jmax,class(GP.rs));


%% 1a. Perform the transform with dmatrix

if exist('dmatrix','var')
    for j = jmax-1:-1:1
        regioncount = nnz(GP.rs(:,j))-1;
        for r = 1:regioncount
            % the index that marks the start of the first subregion
            rs1 = GP.rs(r,j);

            % the index that is one after the end of the second subregion
            rs3 = GP.rs(r+1,j);

            % the number of points in the current region
            n = double(rs3 - rs1);

            if n == 1
                %%% SCALING COEFFICIENT (n == 1)
                dmatrix(rs1,j,:) = dmatrix(rs1,j+1,:);
            elseif n > 1
                % the index that marks the start of the second subregion
                rs2 = rs1+1;
                while rs2 < rs3 && tag(rs2, j+1) ~= 0%%% && rs2 < N+1
                    rs2 = rs2+1;
                end

                % the parent region is a copy of the subregion
                if rs2 == rs3
                    dmatrix(rs1:rs3-1,j,:) = dmatrix(rs1:rs3-1,j+1,:);
                    tag(rs1:rs3-1,j) = tag(rs1:rs3-1,j+1);
                    
                % the parent region has 2 child regions
                else
                    % the number of points in the first subregion
                    n1 = double(rs2 - rs1);

                    % the number of points in the second subregion
                    n2 = double(rs3 - rs2);

                    %%% SCALING COEFFICIENT (n > 1)
                    dmatrix(rs1,j,:)   = ( sqrt(n1)*dmatrix(rs1,j+1,:) + sqrt(n2)*dmatrix(rs2,j+1,:) )/sqrt(n);
                    compinfo(rs1,j) = n1;

                    %%% HAAR-LIKE COEFFICIENT
                    dmatrix(rs1+1,j,:) = ( sqrt(n2)*dmatrix(rs1,j+1,:) - sqrt(n1)*dmatrix(rs2,j+1,:) )/sqrt(n);
                    tag(rs1+1,j) = 1;
                    compinfo(rs1+1,j) = n2;

                    %%% WALSH-LIKE COEFFICIENTS
                    % sweep through the coefficients in subregion 1 and subregion 2

                    % the index of the new coefficient(s) to be created on level j
                    parent = rs1+2;

                    % the index of the current coefficients in subregions 1 and 2
                    child1 = rs1+1;
                    child2 = rs2+1;

                    while parent < rs3
                        % no matching coefficient (use subregion 1)
                        if child1 < rs2 && (child2 == rs3 || tag(child1,j+1) < tag(child2,j+1))
                            dmatrix(parent,j,:) = dmatrix(child1,j+1,:);
                            tag(parent,j) = 2*tag(child1,j+1);
                            child1 = child1+1;
                            parent = parent+1;

                        % no matching coefficient (use subregion 2)
                        elseif child2 < rs3 && (child1 == rs2 || tag(child2,j+1) < tag(child1,j+1))
                            dmatrix(parent,j,:) = dmatrix(child2,j+1,:);
                            tag(parent,j) = 2*tag(child2,j+1);
                            child2 = child2+1;
                            parent = parent+1;

                        % matching coefficients
                        else
                            dmatrix(parent,j,:)   = ( dmatrix(child1,j+1,:) + dmatrix(child2,j+1,:) )/sqrt(2);
                            tag(parent,j)         = 2*tag(child1,j+1);
                            compinfo(parent,j)    = 1;

                            dmatrix(parent+1,j,:) = ( dmatrix(child1,j+1,:) - dmatrix(child2,j+1,:) )/sqrt(2);
                            tag(parent+1,j)       = 2*tag(child1,j+1)+1;
                            compinfo(parent+1,j)  = 1;

                            child1 = child1+1;
                            child2 = child2+1;
                            parent = parent+2;
                        end
                    end
                end
            end
        end
    end
    
%% 1b. Perform the transform without dmatrix

else
    for j = jmax-1:-1:1
        regioncount = nnz(GP.rs(:,j))-1;
        for r = 1:regioncount
            % the index that marks the start of the first subregion
            rs1 = GP.rs(r,j);

            % the index that is one after the end of the second subregion
            rs3 = GP.rs(r+1,j);

            % the number of points in the current region
            n = double(rs3 - rs1);

            if n > 1
                % the index that marks the start of the second subregion
                rs2 = rs1+1;
                while rs2 < rs3 && tag(rs2, j+1) ~= 0
                    rs2 = rs2+1;
                end
                
                % the parent region is a copy of the subregion
                if rs2 == rs3
                    tag(rs1:rs3-1,j) = tag(rs1:rs3-1,j+1);
                    compinfo(rs1:rs3-1,j) = compinfo(rs1:rs3-1,j+1);
                    
                % the parent region has 2 child regions
                else
                    % the number of points in the first subregion
                    n1 = double(rs2 - rs1);

                    % the number of points in the second subregion
                    n2 = double(rs3 - rs2);

                    %%% SCALING COEFFICIENT (n > 1)
                    compinfo(rs1,j) = n1;

                    %%% HAAR-LIKE COEFFICIENT
                    tag(rs1+1,j) = 1;
                    compinfo(rs1+1,j) = n2;

                    %%% WALSH-LIKE COEFFICIENTS
                    % sweep through the coefficients in subregion 1 and subregion 2

                    % the index of the new coefficient(s) to be created on level j
                    parent = rs1+2;

                    % the index of the current coefficients in subregions 1 and 2
                    child1 = rs1+1;
                    child2 = rs2+1;

                    while parent < rs3
                        % no matching coefficient (use subregion 1)
                        if child1 < rs2 && (child2 == rs3 || tag(child1,j+1) < tag(child2,j+1))
                            tag(parent,j) = 2*tag(child1,j+1);
                            child1 = child1+1;
                            parent = parent+1;

                        % no matching coefficient (use subregion 2)
                        elseif child2 < rs3 && (child1 == rs2 || tag(child2,j+1) < tag(child1,j+1))
                            tag(parent,j) = 2*tag(child2,j+1);
                            child2 = child2+1;
                            parent = parent+1;

                        % matching coefficients
                        else
                            tag(parent,j) = 2*tag(child1,j+1);
                            compinfo(parent,j) = 1;

                            tag(parent+1,j) = 2*tag(child1,j+1)+1;
                            compinfo(parent+1,j) = 1;

                            child1 = child1+1;
                            child2 = child2+1;
                            parent = parent+2;
                        end
                    end
                end
            end
        end
    end
end


GP.tag = tag;
GP.compinfo = compinfo;


end