% Make directory of all results files
FileList  = dir(fullfile(pwd,'optResults_*_*.txt'));

% Import and store data from results files
dataStore = [];

efforts = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0]; minDs = [7 15 30 60];
% efforts = [0.5 1.0 2.0 4.0]; minDs = [7 15 30 60];

for j = 1:length(minDs)
    for i = 1:length(efforts)

    filename = sprintf('optResults_%.1f_%d.txt',efforts(i),minDs(j));
    data = importdata(filename);
    dataStore = [dataStore;data]; 
    end
end

% % Read in parameter table
% params = readtable('params.txt');
% P = params.P;
% effort = params.effort;
% minD = params.minD;

%%
% Define colours for IBM colourblind safe palette
IBMblue = [100,143,255]/256;
IBMpurple = [120,93,241]/256;
IBMpink = [221,37,128]/256;
IBMorange = [254,97,0]/256;
IBMyellow = [255,176,0]/256;

%%
close all;

figure(); hold on; box off; set(gca,'LineWidth',2,'FontSize',16);


% Plot month patch stripes
daysPerMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
cD = [0 cumsum(daysPerMonth)];
mids = [16, 14.5, 16, 15.5, 16, 15.5, 16, 16, 15.5, 16, 15.5, 16];
monthMidPoints = cD(1:12) + mids;

x1 = 0.51*ones(1,6); x2 = 4*ones(1,6);
y1 = [cD(2) cD(4) cD(6) cD(8) cD(10) cD(12)]; y2 = [cD(3) cD(5) cD(7) cD(9) cD(11) cD(13)];

p1 = [x1; x2; x2; x1];
p2 = [y1; y1; y2; y2];

mypatches = patch(p1,p2,0.95*[1 1 1]);
mypatches.EdgeAlpha = 0;


% Plot background colour patch regions
patch1 = patch([dataStore(1:8,1);flip(dataStore(1:8,1))],[dataStore(1:8,4);flip(dataStore(1:8,5))],IBMpink);
patch1.LineStyle = 'none'; patch1.FaceAlpha = 0.2;

patch2 = patch([dataStore(9:16,1);flip(dataStore(9:16,1))],[dataStore(9:16,4);flip(dataStore(9:16,5))],IBMorange);
patch2.LineStyle = 'none'; patch2.FaceAlpha = 0.2;

patch3 = patch([dataStore(17:24,1);flip(dataStore(17:24,1))],[dataStore(17:24,4);flip(dataStore(17:24,5))],IBMyellow);
patch3.LineStyle = 'none'; patch3.FaceAlpha = 0.2;

patch4 = patch([dataStore(25:32,1);flip(dataStore(25:32,1))],[dataStore(25:32,4);flip(dataStore(25:32,5))],IBMblue);
patch4.LineStyle = 'none'; patch4.FaceAlpha = 0.2;


% Plot lines
p7s = plot(dataStore(1:8,1),dataStore(1:8,4),'Color',IBMpink,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p7e = plot(dataStore(1:8,1),dataStore(1:8,5),'Color',IBMpink,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

p15s = plot(dataStore(9:16,1),dataStore(9:16,4),'Color',IBMorange,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p15e = plot(dataStore(9:16,1),dataStore(9:16,5),'Color',IBMorange,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

p30s = plot(dataStore(17:24,1),dataStore(17:24,4),'Color',IBMyellow,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p30e = plot(dataStore(17:24,1),dataStore(17:24,5),'Color',IBMyellow,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

p60s = plot(dataStore(25:32,1),dataStore(25:32,4),'Color',IBMblue,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p60e = plot(dataStore(25:32,1),dataStore(25:32,5),'Color',IBMblue,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

ylim([0 365]);
xlabel('Sampling effort (E)');
ylabel('Optimal sampling window [D_{start}^*, D_{end}^*]');

yticks = monthMidPoints;
yticklabels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',''};
set(gca,'YTick',yticks,'YTickLabel',yticklabels,'TickDir','out');




leg = legend([p7s p15s p30s p60s],'\Delta = 7 days','\Delta = 15 days','\Delta = 30 days','\Delta = 60 days','NumColumns',2,'Orientation','vertical');

%%
figure(); hold on; box off; set(gca,'LineWidth',2,'FontSize',16,'YGrid','on');

p7c = plot(dataStore(1:8,1),dataStore(1:8,3),'Color',IBMpink,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p15c = plot(dataStore(9:16,1),dataStore(9:16,3),'Color',IBMorange,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p30c = plot(dataStore(17:24,1),dataStore(17:24,3),'Color',IBMyellow,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p60c = plot(dataStore(25:32,1),dataStore(25:32,3),'Color',IBMblue,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

xlabel('Sampling effort (E)');
ylabel('Change in EDP compared to baseline (%)');

leg = legend([p7c p15c p30c p60c],'\Delta = 7 days','\Delta = 15 days','\Delta = 30 days','\Delta = 60 days','NumColumns',2,'Orientation','vertical');


%% -----------------------------------------------------------------------------------------------
% EXTRACT DATA

% Extract sample size and sample interval vectors
startDayVec = unique(dataStore(:,3))';
endDayVec = unique(dataStore(:,4))';

startSize = length(startDayVec);
endSize = length(endDayVec);

% Extract matrix of the baseline EDP
baselineEDP = reshape(dataStore(:,5),[endSize startSize]);
baselineEDP_smooth = imgaussfilt(baselineEDP,1);

% Extract matrix of the new EDP
newEDP = reshape(dataStore(:,6),[endSize startSize]);
newEDP_smooth = imgaussfilt(newEDP,1);

percChange = 100*(newEDP_smooth-baselineEDP_smooth)./baselineEDP_smooth;

optPoint = min(min(percChange));
[endMin,startMin]=find(percChange==optPoint);
optStart = startDayVec(startMin);
optEnd = endDayVec(endMin);