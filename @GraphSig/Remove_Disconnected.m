function G = Remove_Disconnected(G)
% Remove the nodes from G which are not connected
%
% Input
%   G           a GraphSig object
%
% Output
%   G           a GraphSig object with disconnected nodes removed
%
% 
% 
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



for row = length(G.W):-1:1
    if sum(G.W(row,:)) < 10*eps
        G.W(row,:) = [];
        G.W(:,row) = [];
        if ~isempty(G.xy)
            G.xy(row,:) = [];
        end
        if ~isempty(G.f)
            G.f(row,:) = [];
        end
    end
end


end