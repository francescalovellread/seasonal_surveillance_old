C1 = 0.007265;
a = 0.0004574;
g0 = 0.064913;
tem = 80;

numSims = 100;

day_vals = 1:1:365;
Pvec = C1*exp(g0*(day_vals-tem)-(a/2)*(day_vals-tem).^2);
omega = Pvec/mean(Pvec);

pdf = omega/sum(omega);
cdf = cumtrapz(pdf);
r = rand(1,numSims);

startDays = floor(interp1(cdf,day_vals,r)+1);

temp = [9 9 11 14 19 23 26 26 22 18 13 10];
temp(temp<17)=0;
temp = temp-min(temp); temp = temp/mean(temp);
temp = interp1(1:1:12,temp,linspace(1,12,365));

omega1 = omega;
omega2 = temp;




precip = [46.7 49 45.3 33.9 22.0 13.1 10.8 16.5 42.7 61.6 70.1 58.0];
 
precip2 = [precip(end) precip precip(1)];

precip=precip-min(precip);
    precip=precip/max(precip);
    dryness = 1-precip;
    dryness = [dryness dryness(end)];



% dryness = interp1(1:1:12,dryness,linspace(1,12,365));

daysPerMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
cumulativeDays = [0 cumsum(daysPerMonth)];
mids = [16, 14.5, 16, 15.5, 16, 15.5, 16, 16, 15.5, 16, 15.5, 16];
monthMidPoints = cumulativeDays(1:12) + mids;
cD = cumulativeDays;


daysPerMonth2 = [31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31];
cumulativeDays2 = [0 cumsum(daysPerMonth2)] - 31;
mids2 = [16, 16, 14.5, 16, 15.5, 16, 15.5, 16, 16, 15.5, 16, 15.5, 16, 16];
monthMidPoints2 = cumulativeDays2(1:14) + mids2;


x1 = [cD(2) cD(4) cD(6) cD(8) cD(10) cD(12)]; x2 = [cD(3) cD(5) cD(7) cD(9) cD(11) cD(13)];
y1 = [0 0 0 0 0 0]; y1 = 0.005*ones(1,6);
y2 = [3.5 3.5 3.5 3.5 3.5 3.5];

p1 = [x1; x2; x2; x1];
p2 = [y1; y1; y2; y2];

%%
close all;

figure(); hold on; box off; grid off; set(gca,'LineWidth',2,'FontSize',16,'TickDir','out');

mypatches = patch(p1,p2,0.95*[1 1 1]);
mypatches.EdgeAlpha = 0;

oplot = plot(omega);
oplot.LineWidth = 3;
oplot.Color = [0 0.7 0.6];
meanline = yline(1); meanline.LineStyle = '--'; meanline.LineWidth = 2;
xlim([0 365]);

xlabel('Month'); ylabel('Transmission rate scaling factor: \omega(t)');

% xticks = sort([monthMidPoints cumulativeDays]);
xticks = monthMidPoints;
% xticklabels = {'','Jan','','Feb','','Mar','','Apr','','May','','Jun','','Jul','','Aug','','Sep','','Oct','','Nov','','Dec',''};
xticklabels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

set(gca,'XTick',xticks,'XTickLabel',xticklabels,'XTickLabelRotation',0);
% set(gca,'XTick',xticks,'XTickLabelRotation',0);



%%
figure(); hold on; box off; grid off; set(gca,'LineWidth',2,'FontSize',16,'TickDir','out');

mypatches = patch(p1,p2,0.95*[1 1 1]);
mypatches.EdgeAlpha = 0;

tplot = plot([0 monthMidPoints 365],[0 temp 0]);
tplot.LineWidth = 3;
tplot.Color = [0.5 0.9 0.4];
meanline = yline(1); meanline.LineStyle = '--'; meanline.LineWidth = 2;
xlim([0 365]); ylim([0 2.5]);

xticks = monthMidPoints;
xticklabels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'XTickLabelRotation',0);

xlabel('Month'); ylabel('Symptom onset rate scaling factor: \psi(t)');

%%
figure(); hold on; box off; grid off; set(gca,'LineWidth',2,'FontSize',16,'TickDir','out');

mypatches = patch(p1,p2,0.95*[1 1 1]);
mypatches.EdgeAlpha = 0;

dplot = stairs(cumulativeDays,dryness);
dplot.LineWidth = 3;
dplot.Color = [1 0.6 0];
xlim([0 365]); ylim([0 1]);

xticks = monthMidPoints;
xticklabels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'XTickLabelRotation',0);

xlabel('Month'); ylabel('Seasonal detection sensitivity: \chi(t)');

%%
close all
figure(); hold on; box off; grid off; set(gca,'LineWidth',2,'FontSize',16,'TickDir','out','SortMethod', 'depth');

mypatches = patch(p1,p2,0.95*[1 1 1]);
mypatches.EdgeAlpha = 0;

% mycolour = uisetcolor;
mycolour = [0.3765 0.7451 0.8784];

yyaxis right
pplot = plot(monthMidPoints2,precip2);
pplot.LineStyle = '--';
pplot.LineWidth = 2;
pplot.Marker = 'o'; 
pplot.MarkerFaceColor = mycolour;
pplot.Color = mycolour;
% ylim([0 1]);
ylabel('Average monthly rainfall (mm)');

yyaxis left
dplot = stairs(cumulativeDays,dryness);
dplot.LineWidth = 3;
dplot.Color = [1 0.6 0];
ylim([0 1.001]);
ylabel('Seasonal detection sensitivity: \chi(t)');


xticks = monthMidPoints;
xticklabels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'XTickLabelRotation',0);

xlabel('Month'); xlim([0 365]); 

ax = gca;
% ax.YAxis(1).Color = [1 0.6 0];
ax.YAxis(2).Color = mycolour;
