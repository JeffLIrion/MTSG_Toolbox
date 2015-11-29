% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



fprintf('\n\n');
fprintf('                  N    jmax    HGLET    GHWT coarse-to-fine    GHWT fine-to-coarse\n');


%% Figure 4.1

W = [ 0, 1, 0, 0, 0, 0; ...
      1, 0, 0, 1, 0, 0; ...
      0, 0, 0, 0, 1, 0; ...
      0, 1, 0, 0, 1, 1; ...
      0, 0, 1, 1, 0, 1; ...
      0, 0, 0, 1, 1, 0];

G = GraphSig(W);

GP = PartitionTreeFiedler(G);
[~,rs] = ExtractData(GP);
[~,jmax] = size(rs);
[c2f,f2c] = HGLET_GHWT_BasisCount(GP);
fprintf('Fig 4.1            %d     %d       %d              %d                     %d\n',length(G),jmax-1,c2f,c2f,f2c);


%% Figures 5.1 and 5.2

G = Gpath(6);
W = ExtractData(G);
W(2,3) = 0.1; W(3,2) = 0.1;
G = GraphSig(W);

GP = PartitionTreeFiedler(G);
[~,rs] = ExtractData(GP);
[~,jmax] = size(rs);
[c2f,f2c] = HGLET_GHWT_BasisCount(GP);
fprintf('Figs 5.1 & 5.2     %d     %d       %d              %d                     %d\n',length(G),jmax-1,c2f,c2f,f2c);


%% Fig 3.2

load('Dendrite.mat','G');

% the number of vertices
N = length(G);

% partition the graph
GP = PartitionTreeFiedler(G);

psi = GHWT_BasisVector(G,GP,2,1);
[W,~,f] = ExtractData(psi);
for node = N:-1:1
    if f(node) > 0
        W(node,:) = [];
        W(:,node) = [];
        f(node) = [];
    end
end
G = GraphSig(W);

% partition the new graph (which is a subset of the dendritic tree)
GP = PartitionTreeFiedler(G);

[~,rs] = ExtractData(GP);
[~,jmax] = size(rs);
[c2f,f2c] = HGLET_GHWT_BasisCount(GP);
fprintf('Fig 3.2          %d    %d    %.1e        %.1e              %.1e\n',length(G),jmax-1,c2f,c2f,f2c);


%% Toronto

load('Toronto.mat');

% partition the new graph (which is a subset of the dendritic tree)
GP = PartitionTreeFiedler(G,1);

[~,rs] = ExtractData(GP);
[~,jmax] = size(rs);
[c2f,f2c] = HGLET_GHWT_BasisCount(GP);
fprintf('Toronto         %d    %d     ~10^%d         ~10^%d               ~10^%d\n',length(G),jmax-1,floor(-c2f),floor(-c2f),floor(-f2c));


%% MN road network

load('MN_MutGauss.mat');

% partition the new graph (which is a subset of the dendritic tree)
GP = PartitionTreeFiedler(G);

[~,rs] = ExtractData(GP);
[~,jmax] = size(rs);
[c2f,f2c] = HGLET_GHWT_BasisCount(GP);
fprintf('Minnesota       %d    %d     ~10^%d         ~10^%d               ~10^%d\n\n\n',length(G),jmax-1,floor(-c2f),floor(-c2f),floor(-f2c));


fprintf('\n\nResults may vary slightly due to machine arithmetic\n\n');