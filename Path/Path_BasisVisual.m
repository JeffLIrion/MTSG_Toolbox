function fig = Path_BasisVisual(G,GP,BS,trans,Nmin)
% Display an HGLET/GHWT basis for a path
%
% Input
%   G           a GraphSig object
%   GP          a GraphPart object
%   BS          a BasisSpec object
%   trans       specifies which transform was used for that portion of the
%               signal: 
%                   00 = HGLET with L    ==> blue
%                   01 = HGLET with Lrw  ==> red
%                   10 = HGLET with Lsym ==> black
%                   11 = GHWT            ==> green
%   Nmin        combine segments shorter in length than Nmin
%
% Output
%   fig         the path basis visualization
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% specify trans if necessary
if ~exist('trans','var') || ~islogical(trans)
    levlist = ExtractData(BS);
    trans = false(length(levlist),2);
end

% remove small segments of length less than Nmin
if exist('Nmin','var')
    [GP,BS,trans] = Remove_Small_Segments(GP,BS,trans,Nmin);
end

% extract data
[~,xy,f] = ExtractData(G);
[ind,~] = ExtractData(GP);
[~,cols] = size(xy);
if cols > 1
    fprintf('\n\nThis is not a 1D signal.  Exiting now.\n\n');
    return
end

% extract levlist, levlengths, and coarse-to-fine info
[levlist,levlengths,c2f] = ExtractData(BS,GP);

if ~c2f
    fig = figure;
    BasisVisual(G,GP,BS);
    return
end

% initial plot stuff
fig = figure;
plot(xy(1),f(1),'-b');
hold on
plot(xy(1),f(1),'-r');
plot(xy(1),f(1),'-k');
plot(xy(1),f(1),'-g');
plot(xy,f,'-w');
ylims = ylim;

% the linewidth
lw = 3;

% an index for the x-coordinates
n = 1;

% plot the segments
for row = 1:length(levlist)
    % the indices for the current interval
    if n == 1
        % don't backtrack the first index
        IX = n:n+levlengths(row)-1;
    else
        % backtrack the first index
        IX = n-1:n+levlengths(row)-1;
    end
    
    % HGLET - L
    if ~trans(row,:)
        plot(xy(IX),f(ind(IX)),'-b','LineWidth',lw);
        
    % HGLET - Lrw
    elseif ~trans(row,1) && trans(row,2)
        plot(xy(IX),f(ind(IX)),'-r','LineWidth',lw);
        
    % HGLET - Lsym
    elseif trans(row,1) && ~trans(row,2)
        plot(xy(IX),f(ind(IX)),'-k','LineWidth',lw);
        
    % GHWT
    elseif trans(row,:)
        plot(xy(IX),f(ind(IX)),'-g','LineWidth',lw);
    end
    
    % draw a vertical line
    if row ~= length(levlist)
        x = xy(ind(n+levlengths(row)-1));
        plot([x,x],2*ylims,'-m');
    end
    
    n = n+levlengths(row);    
    
end

xmin = min(xy)-1;
xmax = max(xy)+1;
xlim([xmin xmax])
ylim(ylims)

set(gcf,'color','w');


end




function [GP,BS,trans] = Remove_Small_Segments(GP,BS,trans,Nmin)
% extract data and get constants
[ind,rs] = ExtractData(GP);
[levlist,levlengths] = ExtractData(BS);
levlengths = double(levlengths);
N = length(ind);


% Combine small regions (proceed from smallest on up)
[minlength,row] = min(levlengths);
while minlength < Nmin
    % CASE 1: the segment of minimum length is the first segment
    if row == 1
        % the start and end of the small segment
        rs1 = 1;
        rs2 = rs1+levlengths(row);
        
        % the end of the segment after
        rs3 = rs2+levlengths(row+1);
        
        % absorb the small segment into its neighbor (segment after)
        rs( rs1 < rs(:,levlist(row+1)) & rs(:,levlist(row+1)) < rs3, levlist(row+1)) = rs1;
        levlengths(row+1) = levlengths(row+1)+levlengths(row);
        levlengths(row) = [];
        levlist(row) = [];
        trans(row,:) = [];
        
    % CASE 2: the segment of minimum length is the last segment
    elseif row == length(levlengths)
        % the start and end of the small segment
        rs1 = N+1-levlengths(row);
        rs2 = N+1;
        
        % the start of the segment before
        rs0 = rs1-levlengths(row-1);
        
        % absorb the small segment into its neighbor (segment before)
        rs( rs0 < rs(:,levlist(row-1)) & rs(:,levlist(row-1)) < rs2, levlist(row-1)) = rs2;
        levlengths(row-1) = levlengths(row-1)+levlengths(row);
        levlengths(row) = [];
        levlist(row) = [];
        trans(row,:) = [];
        
    % CASE 3: the segment of minimum length isn't the first or last segment
    else
        % the length of the minimum segment
        n = levlengths(row);
        
        % the start and end of the small segment and its inverse transform
        rs1 = 1+sum(levlengths(1:row-1));
        rs2 = rs1+n;
        
        % the start of the segment before and its inverse transform
        rs0 = rs1-levlengths(row-1);
        
        % the end of the segment after and its inverse transform
        rs3 = rs2+levlengths(row+1);
        
        % the number of its nodes that will be incorporated into the
        % segments to its right and to its left
        add_left = ceil(n/2);
        add_right = n-add_left;
        
        % the new regionstart that is between rs1 and rs2
        rs12 = rs1+add_left;
                
        % absorb the small segment into its neighbors
        rs( rs0  < rs(:,levlist(row-1)) & rs(:,levlist(row-1)) < rs12, levlist(row-1) ) = rs12;
        rs( rs12 < rs(:,levlist(row+1)) & rs(:,levlist(row+1)) < rs3,  levlist(row+1) ) = rs12;
        
        levlengths(row-1) = levlengths(row-1)+add_left;
        levlengths(row+1) = levlengths(row+1)+add_right;
        levlengths(row) = [];
        levlist(row) = [];
        trans(row,:) = [];
    end
    
    % find the new shortest segment
    [minlength,row] = min(levlengths);
    
end

% the new GP and BS
GP = GraphPart(ind,rs);
BS = BasisSpec(levlist);
BS = levlist2levlengths(GP,BS);
end