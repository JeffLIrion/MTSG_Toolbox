function [Vis,GVis] = HGLET_BasisVisual(G,GP,BS,dvec)
% Display an HGLET basis
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



% extract data
[ind,rs] = ExtractData(GP);
[~,jmax] = size(rs);

% extract levlist and levlengths info
[levlist,levlengths] = ExtractData(BS);
if isempty(levlengths)
    BS = levlist2levlengths(GP,BS);
    [levlist,levlengths] = ExtractData(BS);
end

% allocate spece
Vis = zeros(jmax,G.length);


% fill in the matrix
indr0 = 1;
for row = 1:length(levlist)
    indr = indr0:indr0+levlengths(row)-1;
    if exist('dvec','var') == 1
        Vis(levlist(row),indr) = abs(dvec(indr));
    else
        Vis(levlist(row),indr) = 1 + log2(double((1:levlengths(row))'));
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


end