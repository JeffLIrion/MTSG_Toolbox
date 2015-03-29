function G = EditPlotSpecs(G,plotspecs,~)
% Edit the plotspecs of a GraphSig object
%
% Input
%   G           a GraphSig object
%   plotspecs   the new plotspecs
%   ~           if a 3rd argument is given, append the new plotspecs to the
%               current plotspecs
%
% Output
%   G           the GraphSig object with the new plotspecs
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if nargin == 3
    G.plotspecs = strcat(G.plotspecs, ' ', plotspecs);
else
    G.plotspecs = plotspecs;
end


end