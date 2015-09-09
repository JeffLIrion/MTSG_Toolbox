function display(G)
% Display info about the GraphSig ojbect
%
% Input
%   G       a GraphSig object
%
%
% 
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% Input name, Name, and Length
fprintf('\n\nGraphSig object ''%s''\n\n\n',inputname(1));

if ~isempty(G.name)
    fprintf('Name   = %s\n',G.name);
end
if ~isempty(G.plotspecs)
    fprintf('Plot\nSpecs  = %s\n',G.plotspecs);
end
N = length(G.W);
M = nnz(G.W)/2;
fprintf('Nodes  = %d',N);
fprintf('\nEdges  = %d',M);


%% Size and Dim
fprintf('\nDim    = %d',G.dim);


%% [Data]
if issparse(G.W) == 1
    fprintf('\nW      = [%d, %d] sparse %s',G.length,G.length,class(G.W));
else
    fprintf('\nW      = [%d, %d] %s',G.length,G.length,class(G.W));
end

[frows,fcols] = size(G.f);
if issparse(G.f) == 1
    fprintf('\nf      = [%d, %d] sparse %s',frows,fcols,class(G.f));
else
    fprintf('\nf      = [%d, %d] %s',frows,fcols,class(G.f));
end


fprintf('\n\n\n');


%% Plot the data
set(0, 'DefaultFigureVisible', 'on');
GraphSig_Plot(G);
   

end