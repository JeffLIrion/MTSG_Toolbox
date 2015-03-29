function [GP,dmatrixf2c,IX] = FineToCoarse(GP,dmatrix)
% Fill in the fine-to-coarse info (rs2f2c, tagf2c, and compinfof2c) in a
% GraphPart object.  Also, rearrange a matrix of expansion coefficients.  
%
% Input
%   GP          a GraphPart object without fine-to-coarse info (rsf2c,
%               tagf2c, and compinfof2c)
%   dmatrix     a matrix of expansion coefficients in coarse-to-fine
%               arrangement
%
% Output
%   GP          a GraphPart object with fine-to-coarse info
%   dmatrixf2c  a matrix of expansion coefficients in fine-to-coarse
%               arrangement
%   IX          the reordering index for all levels
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminaries

% get constants
[N,jmax] = size(GP.rs);
N = N-1;

% allocate space
GP.rsf2c       = zeros(N+1,jmax,class(GP.rs));
GP.rsf2c(1,:)  = 1;
GP.rsf2c(2,1)  = N+1;
GP.tagf2c      = zeros(N,jmax,class(GP.tag));
GP.compinfof2c = zeros(N,jmax,class(GP.compinfo));
IX = zeros(N,jmax,class(GP.ind));

if nargout == 2
    dmatrixf2c = zeros(size(dmatrix),class(dmatrix));
end

% make sure the coarse-to-fine tag field is filled in
if isempty(GP.tag) || isempty(GP.compinfo)
    GP = GHWT_Core(GP);
end


%% 1. Generate the fine-to-coarse dictionary

% put the coefficients into fine-to-coarse order
for j = 1:jmax
    % put the basis into tag order
    [GP.tagf2c(:,j),IX(:,j)] = sort(GP.tag(:,jmax+1-j));
    GP.compinfof2c(:,j) = GP.compinfo(IX(:,j),jmax+1-j);

    if nargout > 1
        dmatrixf2c(:,j,:) = dmatrix(IX(:,j),jmax+1-j,:);
    end
    
    % fill in the fine-to-coarse regionstarts
    if j > 1
        rsf2crow = 2;
        for row = 2:N
            if GP.tagf2c(row,j) > GP.tagf2c(row-1,j)
                GP.rsf2c(rsf2crow,j) = row;
                rsf2crow = rsf2crow+1;
            end
        end
        GP.rsf2c(rsf2crow,j) = N+1;
    end
end


end