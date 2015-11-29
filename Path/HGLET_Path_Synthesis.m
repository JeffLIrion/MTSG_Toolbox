function [f,GS] = HGLET_Path_Synthesis(dvec,GP,BS,G,~,~)
% Given a vector of HGLET expansion coefficients and info about the graph 
% partitioning and the choice of basis, reconstruct the signal for an
% unweighted path graph
%
% DCT Info:
%  http://vadkudr.org/Algorithms/DTT/DCTI_DCTII/DCTI_DCTII.html#5
%  "Adapted Wavelet Analysis: From Theory to Software" by M.V. Wickerhauser
%  (page 88)
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
[ind,rs] = ExtractData(GP);
[~,jmax] = size(rs);

% fill in the appropriate entries of dmatrix
dmatrix = dvec2dmatrix(dvec,GP,BS);
f = squeeze(dmatrix(:,jmax,:));


%% 1a. Perform the synthesis transform ==> eigenvectors of L

if nargin < 5
    for j = jmax:-1:1
        regioncount = nnz(rs(:,j))-1;
        for r = 1:regioncount
            % the index that marks the start of the region
            rs1 = rs(r,j);

            % the index that is one after the end of the region
            rs3 = rs(r+1,j);

            % the number of points in the current region
            n = double(rs3 - rs1);

            % only proceed forward if coefficients do not exist
            if (j == jmax || nnz(dmatrix(rs1:rs3-1,j+1,:)) == 0) && nnz(dmatrix(rs1:rs3-1,j,:)) > 0

                if n == 1
                    f(rs1,:) = dmatrix(rs1,j,:);
                elseif n > 1
                    % reconstruct the signal ==> iDCT-II
                    f(rs1:rs3-1,:) = idct(squeeze(dmatrix(rs1:rs3-1,j,:)));
                end

            end
        end
    end


%% 1b. Perform the synthesis transform ==> eigenvectors of L_rw

elseif nargin == 5
    for j = jmax:-1:1
        regioncount = nnz(rs(:,j))-1;
        for r = 1:regioncount
            % the index that marks the start of the region
            rs1 = rs(r,j);

            % the index that is one after the end of the region
            rs3 = rs(r+1,j);

            % the number of points in the current region
            n = double(rs3 - rs1);

            % only proceed forward if coefficients do not exist
            if (j == jmax || nnz(dmatrix(rs1:rs3-1,j+1,:)) == 0)  && nnz(dmatrix(rs1:rs3-1,j,:)) > 0

                if n == 1
                    f(rs1,:) = dmatrix(rs1,j,:);
                elseif n > 1
                    % reconstruct the signal ==> iDCT-I
                    D = 2*ones(n,1);
                    D([1,n]) = 1;
                    f(rs1:rs3-1,:) = bsxfun(@times, D.^(-0.5), idct1(squeeze(dmatrix(rs1:rs3-1,j,:))));
                end

            end
        end
    end


%% 1c. Perform the synthesis transform ==> eigenvectors of L_sym

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

            % only proceed forward if expansion coefficients exist and
            % signal values do not
            if (j == jmax || nnz(dmatrix(rs1:rs3-1,j+1,:)) == 0)  && nnz(dmatrix(rs1:rs3-1,j,:)) > 0

                if n == 1
                    f(rs1,:) = dmatrix(rs1,j,:);
                elseif n > 1
                    % reconstruct the signal ==> iDCT-I
                    f(rs1:rs3-1,:) = idct1(squeeze(dmatrix(rs1:rs3-1,j,:)));
                end
            end
        end
    end
end



% put the reconstructed values in the correct order
f(ind,:) = f;


% create a GraphSig object with the reconstructed data
if nargout == 2
    if exist('G','var') == 1 && isa(G,'GraphSig')
        GS = ReplaceData(G,f);
    else
        GS = 0;
    end
end


end
