function G = ReplaceData(G,newdata)
% Replace the data in a GraphSig object
%
% Input
%   G           the input GraphSig object whose data is to be replaced
%   newdata     the new data (either a vector or a GraphSig object)
%
%
% Output
%   G           the new GraphSig object with the new data vector
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if strcmp(class(newdata),'GraphSig') == 1
    G.f = newdata.f;
else
    [rows,cols] = size(newdata);
    if G.length == rows
        G.f = newdata;
    elseif G.length == cols
        G.f = newdata';
    elseif isscalar(newdata)
        G.f = newdata*ones(G.length,1);
    end
end


end