% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if exist('WavePath.m','file')
    load('Blocks.mat');
    G0 = G;
    [~,~,f0] = ExtractData(G0);
    
    load('Blocks_Noisy.mat');
    [~,~,f] = ExtractData(G);
    N = length(G);
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
    fprintf('\n\nBlocks:\n');
    fprintf('N = %d\n',length(f));
    fprintf('Original SNR = 11.95 dB\n');
    fprintf('Final SNR = %.2f dB\n',SNR(f0,f2));
else
    fprintf('\n\nPlease install WaveLab if you would like to generate this figure.\n\n<a href="http://statweb.stanford.edu/~wavelab/">http://statweb.stanford.edu/~wavelab/</a>\n\n');
end