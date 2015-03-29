% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



N = 6;


%% 0.  Create a GraphSig object for a path graph of length N

% the weight matrix
e = ones(N,1);
W = spdiags([ e 0*e e], [-1 0 1],N,N);

% the spatial coordinates
xy = (1:N)';

% the signal on the graph
f = xy .* sin(xy);

% additional info for the GraphSig object
dataname = sprintf('A signal on the path graph of length %d',N);
plotspecs = 'notitle stem';

% create the GraphSig object
G = GraphSig(W,xy,f,dataname,plotspecs);
clear N e W xy f dataname plotspecs
display(G);
pause



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% HGLET Demo %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.1.  Partition the graph using the unnormalized Laplacian

GP = PartitionTreeFiedler(G);


%% 1.2.  Perform the HGLET transform on the graph data

dmatrix = HGLET_Analysis(G,GP);


%% 1.3.  Select the best-basis and display the chosen basis

p = 0.1;
[dvec,BS] = HGLET_BestBasis(dmatrix,GP,p);

[~,GVis] = HGLET_BasisVisual(G,GP,BS);
GVis = EditName(GVis,'Selected HGLET Best Basis');
GVis = EditPlotSpecs(GVis,'red',1);
GraphSig_Overlay(GVis,G,2,'-k',1,'-c');
pause


%% 1.4.  Reconstruct the signal using the best-basis expansion coefficients

% create a GraphSig object containing only the weight matrix
W = ExtractData(G); GW = GraphSig(W);

% reconstruct the signal
[~,GS] = HGLET_Synthesis(dvec,GP,BS,GW);

% check the relative error to make sure it is reconstructed up to machine precision
fprintf('\n\nHGLET relative reconstruction error = %.2e\n\n',norm(G-GS,'fro')/norm(G,'fro'));
pause
clear GP dmatrix dvec BS p W GW GS



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% GHWT Demo %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2.1.  Partition the graph using the random walk Laplacian

GP0 = PartitionTreeFiedler(G,1);
display(GP0);


%% 2.2.  Perform the GHWT transform on the graph data

[dmatrix,GP] = GHWT_Analysis(G,GP0);
display(GP);
fprintf('\n\nNote that GP contains "tag" and "compinfo" data, whereas GP0 does not.\n\n');
clear GP0
pause


%% 2.3.  Select the best-basis

p = 0.1;
[~,BS] = GHWT_BestBasis(dmatrix,GP,p);
display(BS);
[~,~,c2f] = ExtractData(BS);
if c2f
    fprintf('\n\nNote that the best-basis is from the coarse-to-fine dictionary.\n\n');
else
    fprintf('\n\nNote that the best-basis is from the fine-to-coarse dictionary.\n\n');
end
clear dmatrix
pause


%% 2.4.  Generate a new signal on the graph

[~,xy] = ExtractData(G);
f2 = xy.^2;
G2 = ReplaceData(G,f2);
display(G2);
pause


%% 2.5.  Perform the GHWT transform on the new graph data

dmatrix = GHWT_Analysis(G2,GP);
clear G2 xy


%% 2.6.  Select the expansion coefficients corresponding to the best-basis for G

% select the expansion coefficients
dvec = dmatrix2dvec(dmatrix,GP,BS);


%% 2.7.  Reconstruct the signal from G2

f2S = GHWT_Synthesis(dvec,GP,BS);

% check the relative error to make sure it is reconstructed up to machine precision
fprintf('\n\nGHWT relative reconstruction error = %.2e\n\n',norm(f2-f2S,2)/norm(f2,2));
clear G GP p BS c2f f2 f2S dmatrix dvec