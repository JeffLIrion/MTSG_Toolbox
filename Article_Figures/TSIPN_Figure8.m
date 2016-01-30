% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% specify linewidth and fontsize
lw = 2;
fs = 12;

% load the noisy signal
load('Toronto_SNR7.mat');

% partition the graph
GP = PartitionTreeFiedler(G,1);

% analyze the noisy signal
dH = HGLET_Analysis(G,GP);

% find the best basis for the noisy signal
[dvec,BS,~,tau] = HGLET_GHWT_BestBasis_minrelerror(dH,0,0,0,GP,G);

% the magnitudes of the coefficients
x = sort(abs(dvec),'descend');

% compute the relative error curve
r = orth2relerror(dvec);

% shift the data points
x = [x; 0];
r = [r(1); r];

% reorder the data points
xmax = max(x);
ymax = max(r);
x = x(end:-1:1)/xmax;
y = r(end:-1:1)/ymax;


%% first time

% a unit vector pointing from (x1,y1) to (xN,yN)
v = [x(end)-x(1), y(end)-y(1)];
v = v/norm(v,2);

% subtract (x1,y1) from the coordinates
xy = [x-x(1), y-y(1)];

% the hypotenuse
H = (sum(xy.^2,2)).^0.5;

% the adjacent side
A = xy*v';

% the opposite side
O = (H.^2-A.^2).^0.5;

[~,IX1] = max(O);

close all
plot(x,y,'LineWidth',lw)
axis square;
ylim([0,1.004]);
xlim([0,1.004]);
hold on

% plot the line
plot([x(1),1.004*x(end)],[y(1),1.004*y(end)],'-g','LineWidth',lw)

% drop the perpendicular
dropto = [x(1), y(1)]+A(IX1)*v;
plot([x(IX1),dropto(1)],[y(IX1),dropto(2)],'-g','LineWidth',lw)


%% second time

% a unit vector pointing from (x1,y1) to (x(IX1),y(IX1))
v = [x(IX1)-x(1), y(IX1)-y(1)];
v = v/norm(v,2);

% subtract (x1,y1) from the coordinates
xy = [x(1:IX1)-x(1), y(1:IX1)-y(1)];

% the hypotenuse
H = (sum(xy.^2,2)).^0.5;

% the adjacent side
A = xy*v';

% the opposite side
O = (H.^2-A.^2).^0.5;

[~,IX2] = max(O);

% plot the line
plot([x(1),x(IX1)],[y(1),y(IX1)],'-r','LineWidth',lw)

% drop the perpendicular
dropto = [x(1), y(1)]+A(IX2)*v;
plot([x(IX2),dropto(1)],[y(IX2),dropto(2)],'-r','LineWidth',lw)

% plot the points on the relative error curve
scatter(x(IX1),y(IX1),'go','MarkerFaceColor','g')
scatter(x(IX2),y(IX2),'ro','MarkerFaceColor','r')

% the tick mark labels
set(gca,'XTick',0:0.251:1.004);
set(gca,'YTick',0:0.251:1.004);
XTL = get(gca,'XTickLabel');
YTL = get(gca,'YTickLabel');
if iscell(XTL)
    XTL = str2num(char(XTL));
    YTL = str2num(char(YTL));
else
    XTL = str2num(XTL);
    YTL = str2num(YTL);
end
set(gca,'XTickLabel',XTL*xmax);
set(gca,'YTickLabel',num2str(YTL*ymax,'%.2f'));

set(gcf,'color','w');
xlabel('Threshold','FontSize',12);
ylabel('Relative Error','FontSize',12);