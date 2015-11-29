% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~exist('PieceRegular.mat','file')
    fprintf('\n\n\nIn order to generate these results, install WaveLab (<a href="http://statweb.stanford.edu/~wavelab/">http://statweb.stanford.edu/~wavelab/</a>),\nrun WavePath.m, then run this script again.\n\n\n');
    return
end

if exist('WavePath.m','file')
    load('PieceRegular.mat','G');
    [~,~,f0] = ExtractData(G);
    f0 = [f0(1); f0; f0(end-1:end)];
    G0 = G;

    load('PieceRegular_SNR20.mat','G');
    [~,~,f] = ExtractData(G);
    f = [f(1); f; f(end-1:end)];
    N = length(f);
    x = (1:N)';
    
    QMF_Filter = MakeONFilter('Symmlet',8);
	L = 5;
    thr = sqrt(2* log(N));    
    
    xhblocks = zeros(size(f));
    tiwtblocks = FWT_TI(f,L,QMF_Filter);
    
    [nrow,ncol] = size(tiwtblocks);
    thrwtblocks = tiwtblocks;
    thrwtblocks(:,2:ncol) = SoftThresh(thrwtblocks(:,2:ncol),thr);
    
    f2 = IWT_TI(thrwtblocks,QMF_Filter)';
    
    figure;
    plot(x,f,'-b','LineWidth',3);
    ylims = ylim;
    close(gcf);
    
    figure;
    plot(x,f2,'-b','LineWidth',3);
    ylim(ylims);
    xlim([0,N+1]);
    set(gcf,'color','w');
    
    % info about the run
    fprintf('\n\nPiece-Regular:\n');
    fprintf('N = %d\n',length(f));
    fprintf('Original SNR = 20.00 dB\n');
    fprintf('Final SNR = %.2f dB\n',SNR(f0,f2));
else
    fprintf('\n\nPlease install WaveLab if you would like to generate this figure.\n\n<a href="http://statweb.stanford.edu/~wavelab/">http://statweb.stanford.edu/~wavelab/</a>\n\n');
end