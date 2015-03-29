function [dvec,kept] = dvec_Threshold(dvec,SORH,keep,GP,BS)
% Threshold HGLET / GHWT coefficients
%
% Inputs
%   dvec        the vector of expansion coefficients
%   SORH        use soft ('s') or hard ('h') thresholding
%   keep        either an integer or a fraction between 0 and 1 which says
%               how many coefficients should be kept
%   GP          a GraphPart object, used to identify scaling coefficients
%   BS          a BasisSpec object, used to identify scaling coefficients
%
% Outputs
%   dvec        the thresholded expansion coefficients
%   kept        the number of coefficients kept
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminaries
N = length(dvec);

% if keep is given as a fraction
if keep < 1-eps && keep > eps
    kept = round(keep*N);
    
% if keep is given as an integer
else
    kept = round(keep);
    if kept > N
        kept = N;
    elseif kept < 0
        kept = 1;
    end
end

% if keep is too large
if kept >= N
    kept = N;
    return
end


%% 1. Threshold the coefficients
if strcmpi(SORH,'s') || strcmpi(SORH,'soft')
    [~,ind] = sort(abs(dvec),'descend');
    T = abs(dvec(ind(kept+1)));
    
    % if BS is given, then use tag information so that scaling coefficients
    % are not soft-thresholded
    if nargin == 5
        % extract the tag info
        [~,~,tag] = ExtractData(GP);
        if isempty(tag)
            GP = GHWT_Core(GP);
            [~,~,tag] = ExtractData(GP);
        end
        
        % convert the tag matrix to a vector corresponding to BS
        tag = dmatrix2dvec(tag,GP,BS);
        
        % soft-threshold the non-scaling coefficients
        ind = (1:N)';
        ind = ind( tag > 0 );
        dvec(ind) = sign( dvec(ind) ) .* (abs( dvec(ind) )-T) .* (abs( dvec(ind) )-T > 0);
        
    % if BS is not given, soft-threshold all coefficients
    else
        dvec = sign(dvec) .* (abs(dvec)-T) .* (abs(dvec)-T > 0);
    end
    
    kept = nnz(dvec);

elseif strcmpi(SORH,'h') || strcmpi(SORH,'hard')
    [~,ind] = sort(abs(dvec),'descend');
    dvec(ind(kept+1:N)) = 0;
    kept = nnz(dvec);

else
    kept = N;
end


end