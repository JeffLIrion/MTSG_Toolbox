% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% the data sets and corresponding signal-to-noise ratios
datasets = {'MN_MutGauss','Dendrite','Toronto'};
SNRs = [5;8;7];

% for storing results
tau1 = zeros(3,1);
tau3 = zeros(3,1);
tau5 = zeros(3,1);
results = zeros(6,3);

for t = 1:3
    %% load and process the original and noisy data 
    
    % load the original data set
    load(sprintf('%s.mat',datasets{t}));
    if strcmp(datasets{t},'MN_MutGauss')
        G = Adj2InvEuc(G);
    end
    
    % info about the data set
    N = length(G);
    [~,~,f] = ExtractData(G);
    
    % partition the graph
    GP = PartitionTreeFiedler(G,1);

    % analyze the true signal
    [dH,dHrw,dHsym] = HGLET_Analysis_All(G,GP);
    dG = GHWT_Analysis(G,GP);    

    
    % load the noisy data set
    load(sprintf('%s_SNR%d.mat',datasets{t},SNRs(t)));

    % analyze the noisy signal
    [dmatrixH,dmatrixHrw,dmatrixHsym] = HGLET_Analysis_All(G,GP);
    dmatrixG = GHWT_Analysis(G,GP);

    
    %% find best bases and compute relative error & SNR curves
    
    % 1) HGLET best basis with L
    [dvec1,BS1,trans1,tau1(t)] = HGLET_GHWT_BestBasis_minrelerror(dmatrixH,0,0,0,GP,G);
    dvec1_true = dmatrix2dvec(dH,GP,BS1);
    r1 = orth2relerror(dvec1);
    s1 = orth2SNR(dvec1,dvec1_true,GP,BS1,'s');
    x1 = sort(abs(dvec1),'descend');

    % 2) Laplacian eigenvectors
    BS2 = LevelBasisSpec(GP,0);
    dvec2 = dmatrix2dvec(dmatrixH,GP,BS2);
    dvec2_true = dmatrix2dvec(dH,GP,BS2);
    r2 = orth2relerror(dvec2);
    s2 = orth2SNR(dvec2,dvec2_true,GP,BS2,'s');
    x2 = sort(abs(dvec2),'descend');

    % 3) GHWT best basis
    [dvec3,BS3,trans3,tau3(t)] = HGLET_GHWT_BestBasis_minrelerror(0,0,0,dmatrixG,GP,G);
    dvec3_true = dmatrix2dvec(dG,GP,BS3);
    r3 = orth2relerror(dvec3);
    s3 = orth2SNR(dvec3,dvec3_true,GP,BS3,'s');
    x3 = sort(abs(dvec3),'descend');

    % 4) Haar basis
    BS4 = HaarBasisSpec(GP);
    dvec4 = dmatrix2dvec(dmatrixG,GP,BS4);
    dvec4_true = dmatrix2dvec(dG,GP,BS4);
    r4 = orth2relerror(dvec4);
    s4 = orth2SNR(dvec4,dvec4_true,GP,BS4,'s');
    x4 = sort(abs(dvec4),'descend');

    % 5) Hybrid best basis (HGLET with L/Lrw/Lsym + GHWT coarse-to-fine)
    [dvec5,BS5,trans5,tau5(t)] = HGLET_GHWT_BestBasis_minrelerror(dmatrixH,dmatrixHrw,dmatrixHsym,dmatrixG,GP,G,1);
    B5 = HGLET_GHWT_Synthesis(eye(N),GP,BS5,trans5,G);
    dvec5_true = dmatrices2dvec(dH,dHrw,dHsym,dG,GP,BS5,trans5);
    r5 = nonorth2relerror(dvec5,B5);
    s5 = nonorth2SNR(dvec5,f,GP,BS5,'s',B5);
    x5 = sort(abs(dvec5),'descend');
    
    % 6) Walsh basis (GHWT level j=0)
    BS6 = LevelBasisSpec(GP,0);
    dvec6 = dmatrix2dvec(dmatrixG,GP,BS6);
    dvec6_true = dmatrix2dvec(dG,GP,BS6);
    r6 = orth2relerror(dvec6);
    s6 = orth2SNR(dvec6,dvec6_true,GP,BS6,'s');
    x6 = sort(abs(dvec6),'descend');


    %% modify the coefficient magnitudes, relative errors, and SNRs
    x1 = [x1; 0];
    x2 = [x2; 0];
    x3 = [x3; 0];
    x4 = [x4; 0];
    x5 = [x5; 0];
    x6 = [x6; 0];

    r1 = [r1(1); r1];
    r2 = [r2(1); r2];
    r3 = [r3(1); r3];
    r4 = [r4(1); r4];
    r5 = [r5(1); r5];
    r6 = [r6(1); r6];

    s1 = [s1(1); s1];
    s2 = [s2(1); s2];
    s3 = [s3(1); s3];
    s4 = [s4(1); s4];
    s5 = [s5(1); s5];
    s6 = [s6(1); s6];


    %% denoise via elbow

    idx1 = find_elbow(x1,r1);
    idx2 = find_elbow(x2,r2);
    idx3 = find_elbow(x3,r3);
    idx4 = find_elbow(x4,r4);
    idx5 = find_elbow(x5,r5);
    idx6 = find_elbow(x6,r6);

    idx1 = idx1-1+find_elbow(x1(idx1:end),r1(idx1:end));
    idx2 = idx2-1+find_elbow(x2(idx2:end),r2(idx2:end));
    idx3 = idx3-1+find_elbow(x3(idx3:end),r3(idx3:end));
    idx4 = idx4-1+find_elbow(x4(idx4:end),r4(idx4:end));
    idx5 = idx5-1+find_elbow(x5(idx5:end),r5(idx5:end));
    idx6 = idx6-1+find_elbow(x6(idx6:end),r6(idx6:end));

    dvec1 = dvec_Threshold(dvec1,'s',idx1,GP,BS1);
    dvec2 = dvec_Threshold(dvec2,'s',idx2,GP,BS2);
    dvec3 = dvec_Threshold(dvec3,'s',idx3,GP,BS3);
    dvec4 = dvec_Threshold(dvec4,'s',idx4,GP,BS4);
    dvec5 = dvec_Threshold(dvec5,'s',idx5,GP,BS5);
    dvec6 = dvec_Threshold(dvec6,'s',idx6,GP,BS6);
    
    results(1,t) = SNR(GraphSig(sparse(N,N),[],dvec1_true),GraphSig(sparse(N,N),[],dvec1));
    results(2,t) = SNR(GraphSig(sparse(N,N),[],dvec2_true),GraphSig(sparse(N,N),[],dvec2));
    results(3,t) = SNR(GraphSig(sparse(N,N),[],dvec3_true),GraphSig(sparse(N,N),[],dvec3));
    results(4,t) = SNR(GraphSig(sparse(N,N),[],dvec4_true),GraphSig(sparse(N,N),[],dvec4));
    results(5,t) = SNR(GraphSig(sparse(N,N),[],B5*dvec5_true),GraphSig(sparse(N,N),[],B5*dvec5));
    results(6,t) = SNR(GraphSig(sparse(N,N),[],dvec6_true),GraphSig(sparse(N,N),[],dvec6));

end


% print a table of the results
fprintf('\n\n(Results may differ depending on your machine)\n\n');
fprintf('                           MN_MutGauss          Dendrite             Toronto\n');
fprintf('                            (%.2f dB)           (%.2f dB)           (%.2f dB)\n\n',SNRs(1),SNRs(2),SNRs(3));
fprintf('HGLET (L) Best Basis     %5.2f dB (t=%.1f)    %5.2f dB (t=%.1f)    %5.2f dB (t=%.1f)\n',results(1,1), tau1(1), results(1,2), tau1(2), results(1,3), tau1(3));
fprintf('Lap. Eigenvectors (L)    %5.2f dB            %5.2f dB            %5.2f dB\n',results(2,1), results(2,2), results(2,3));
fprintf('GHWT Best Basis          %5.2f dB (t=%.1f)    %5.2f dB (t=%.1f)    %5.2f dB (t=%.1f)\n',results(3,1), tau3(1), results(3,2), tau3(2), results(3,3), tau3(3));
fprintf('Haar                     %5.2f dB            %5.2f dB            %5.2f dB\n',results(4,1), results(4,2), results(4,3));
fprintf('Hybrid Best Basis        %5.2f dB (t=%.1f)    %5.2f dB (t=%.1f)    %5.2f dB (t=%.1f)\n',results(5,1), tau5(1), results(5,2), tau5(2), results(5,3), tau5(3));
fprintf('Walsh Basis              %5.2f dB            %5.2f dB            %5.2f dB\n',results(6,1), results(6,2), results(6,3));
fprintf('\n\n');