% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% specify linewidth and fontsize
lw = 2;
fs = 12;

% load the noisy signal
load('MN_MutGauss_SNR5.mat');
G = Adj2InvEuc(G);

% partition the graph
GP = PartitionTreeFiedler(G,1);

% analyze the true and noisy signals
dG = GHWT_Analysis(G,GP);

% find the best basis for the noisy signal
[dvec,BS,~,tau] = HGLET_GHWT_BestBasis_minrelerror(0,0,0,dG,GP,G);

% the magnitudes of the coefficients
x = sort(abs(dvec),'descend');

% compute the relative error curve
r = orth2relerror(dvec);

% shift the data points
x = [x; 0];
r = [r(1); r];

% reorder the data points
x = x(end:-1:1)/max(x);
y = r(end:-1:1)/max(r);


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
ylim([0,1]);
xlim([0,max(x)]);
hold on

% plot the line
plot([x(1),x(end)],[y(1),y(end)],'-r','LineWidth',lw)

% drop the perpendicular
dropto = [x(1), y(1)]+A(IX1)*v;
plot([x(IX1),dropto(1)],[y(IX1),dropto(2)],'-r','LineWidth',lw)


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
plot([x(1),x(IX1)],[y(1),y(IX1)],'-g','LineWidth',lw)

% drop the perpendicular
dropto = [x(1), y(1)]+A(IX2)*v;
plot([x(IX2),dropto(1)],[y(IX2),dropto(2)],'-g','LineWidth',lw)
scatter(x(IX2),y(IX2),'go','MarkerFaceColor','g')
scatter(x(IX1),y(IX1),'ro','MarkerFaceColor','r')

set(gcf,'color','w');
xlabel('(Scaled) Threshold','FontSize',12);
ylabel('(Scaled) Relative Error','FontSize',12);