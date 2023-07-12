
%% -----------------------------------------------------------------------------------------------
% READ IN DATA FROM TWO SPECIES WINDOW CASE 

cd '/Users/francescalovell-read/Documents/DPhil/Chapter_4/CHAPTER_4_CODE/Two_species_window/ARC_original/OPTIMAL_RESULTS'

% Make directory of all results files
FileList  = dir(fullfile(pwd,'optResults_*_*.txt'));

% Import and store data from results files
dataStore2S = [];

efforts = [0.5 1.0 2.0 4.0]; minDs = [7.0 15.0 30.0 60.0]; minDs = [7 15 30 60];
efforts = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0]; minDs = [7 15 30 60];

for j = 1:length(minDs)
    for i = 1:length(efforts)

    filename = sprintf('optResults_%.1f_%d.txt',efforts(i),minDs(j));
    data = importdata(filename);
    dataStore2S = [dataStore2S;data]; 
    end
end


%% -----------------------------------------------------------------------------------------------
% READ IN DATA FROM ONE SPECIES WINDOW CASE (FOR COMPARISON LATER)

cd '/Users/francescalovell-read/Documents/DPhil/Chapter_4/CHAPTER_4_CODE/One_species_window/ARC_original/OPTIMAL_RESULTS'

% Make directory of all results files
FileList  = dir(fullfile(pwd,'optResults_*_*.txt'));

% Import and store data from results files
dataStore1S = [];

efforts = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0]; minDs = [7 15 30 60];
% efforts = [0.5 1.0 2.0 4.0]; minDs = [7 15 30 60];

for j = 1:length(minDs)
    for i = 1:length(efforts)

    filename = sprintf('optResults_%.1f_%d.txt',efforts(i),minDs(j));
    data = importdata(filename);
    dataStore1S = [dataStore1S;data]; 
    end
end


%% -----------------------------------------------------------------------------------------------
% Define colours for IBM colourblind safe palette

IBMblue = [100,143,255]/256;
IBMpurple = [120,93,241]/256;
IBMpink = [221,37,128]/256;
IBMorange = [254,97,0]/256;
IBMyellow = [255,176,0]/256;


%% -----------------------------------------------------------------------------------------------
% PLOT FIGURE SHOWING OPTIMAL WINDOWS IN TWO SPECIES CASE

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
patch1 = patch([dataStore2S(1:8,1);flip(dataStore2S(1:8,1))],[dataStore2S(1:8,5);flip(dataStore2S(1:8,6))],IBMpink);
patch1.LineStyle = 'none'; patch1.FaceAlpha = 0.2;

patch2 = patch([dataStore2S(9:16,1);flip(dataStore2S(9:16,1))],[dataStore2S(9:16,5);flip(dataStore2S(9:16,6))],IBMorange);
patch2.LineStyle = 'none'; patch2.FaceAlpha = 0.2;

patch3 = patch([dataStore2S(17:24,1);flip(dataStore2S(17:24,1))],[dataStore2S(17:24,5);flip(dataStore2S(17:24,6))],IBMyellow);
patch3.LineStyle = 'none'; patch3.FaceAlpha = 0.2;

patch4 = patch([dataStore2S(25:32,1);flip(dataStore2S(25:32,1))],[dataStore2S(25:32,5);flip(dataStore2S(25:32,6))],IBMblue);
patch4.LineStyle = 'none'; patch4.FaceAlpha = 0.2;


% Plot lines
p7s = plot(dataStore2S(1:8,1),dataStore2S(1:8,5),'Color',IBMpink,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p7e = plot(dataStore2S(1:8,1),dataStore2S(1:8,6),'Color',IBMpink,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

p15s = plot(dataStore2S(9:16,1),dataStore2S(9:16,5),'Color',IBMorange,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p15e = plot(dataStore2S(9:16,1),dataStore2S(9:16,6),'Color',IBMorange,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

p30s = plot(dataStore2S(17:24,1),dataStore2S(17:24,5),'Color',IBMyellow,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p30e = plot(dataStore2S(17:24,1),dataStore2S(17:24,6),'Color',IBMyellow,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

p60s = plot(dataStore2S(25:32,1),dataStore2S(25:32,5),'Color',IBMblue,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p60e = plot(dataStore2S(25:32,1),dataStore2S(25:32,6),'Color',IBMblue,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);


xlabel('Sampling effort (E)');
ylabel('Optimal sampling window [D_{start}^*, D_{end}^*]');

ylim([0 365]);

yticks = monthMidPoints;
yticklabels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',''};
set(gca,'YTick',yticks,'YTickLabel',yticklabels,'TickDir','out');


leg = legend([p7s p15s p30s p60s],'\Delta = 7 days','\Delta = 15 days','\Delta = 30 days','\Delta = 60 days','NumColumns',2,'Orientation','vertical');
% leg = legend([p7s p15s p30s],'minD = 7 days','minD = 15 days','minD = 30 days','NumColumns',2,'Orientation','vertical');


%% -----------------------------------------------------------------------------------------------
% PLOT FIGURE SHOWING EDP REDUCTIONS COMPARED TO UNIFORM SEASONAL BASELINE

figure(); hold on; box off; set(gca,'LineWidth',2,'FontSize',16,'YGrid','on');

p7c = plot(dataStore2S(1:8,1),dataStore2S(1:8,3),'Color',IBMpink,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p15c = plot(dataStore2S(9:16,1),dataStore2S(9:16,3),'Color',IBMorange,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p30c = plot(dataStore2S(17:24,1),dataStore2S(17:24,3),'Color',IBMyellow,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p60c = plot(dataStore2S(25:32,1),dataStore2S(25:32,3),'Color',IBMblue,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

xlabel('Sampling effort (E)');
ylabel('Change in EDP compared to USB (%)');

leg = legend([p7c p15c p30c p60c],'\Delta = 7 days','\Delta = 15 days','\Delta = 30 days','\Delta = 60 days','NumColumns',2,'Orientation','vertical','Location','southeast');
% leg = legend([p7c p15c p30c],'minD = 7 days','minD = 15 days','minD = 30 days','NumColumns',2,'Orientation','vertical','Location','southeast');


%% -----------------------------------------------------------------------------------------------
% PLOT FIGURE SHOWING EDP REDUCTIONS COMPARED TO WINDOW SEASONAL BASELINE

changes = (dataStore2S(:,3)-dataStore1S(:,3));
% changes = -100*(dataStore2S(:,3)-dataStore1S(:,3))./dataStore1S(:,3);


figure(); hold on; box off; set(gca,'LineWidth',2,'FontSize',16,'YGrid','on');

p7c = plot(dataStore2S(1:8,1),changes(1:8),'Color',IBMpink,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p15c = plot(dataStore2S(9:16,1),changes(9:16),'Color',IBMorange,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p30c = plot(dataStore2S(17:24,1),changes(17:24),'Color',IBMyellow,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p60c = plot(dataStore2S(25:32,1),changes(25:32),'Color',IBMblue,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

xlabel('Sampling effort (E)');
ylabel('Change in EDP compared to WSB (%)');

leg = legend([p7c p15c p30c p60c],'\Delta = 7 days','\Delta = 15 days','\Delta = 30 days','\Delta = 60 days','NumColumns',2,'Orientation','vertical','Location','southeast');
% leg = legend([p7c p15c p30c],'minD = 7 days','minD = 15 days','minD = 30 days','NumColumns',2,'Orientation','vertical','Location','southeast');


%% -----------------------------------------------------------------------------------------------
% PLOT FIGURE SHOWING OPTIMAL NUMBER OF SENTINELS

figure(); hold on; box off; set(gca,'LineWidth',2,'FontSize',16,'YGrid','on');

p7c = plot(dataStore2S(1:8,1),dataStore2S(1:8,7),'Color',IBMpink,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p15c = plot(dataStore2S(9:16,1),dataStore2S(9:16,7),'Color',IBMorange,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p30c = plot(dataStore2S(17:24,1),dataStore2S(17:24,7),'Color',IBMyellow,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p60c = plot(dataStore2S(25:32,1),dataStore2S(25:32,7),'Color',IBMblue,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

xlabel('Sampling effort (E)');
ylabel('Optimal number of sentinels (P_S^*)');

leg = legend([p7c p15c p30c p60c],'\Delta = 7 days','\Delta = 15 days','\Delta = 30 days','\Delta = 60 days','NumColumns',2,'Orientation','vertical','Location','southwest');
% leg = legend([p7c p15c p30c],'minD = 7 days','minD = 15 days','minD = 30 days','NumColumns',2,'Orientation','vertical','Location','southeast');


%%
figure(); hold on; box off; set(gca,'LineWidth',2,'FontSize',16,'YGrid','on');

sampleSizes = ((dataStore2S(1:32,1).*dataStore2S(1:32,2))./((dataStore2S(1:32,6)-dataStore2S(1:32,5)+1)/365));

NSvals = dataStore2S(1:32,8)./min(sampleSizes(1:32),dataStore2S(1:32,7));
NSvals(NSvals>1)=1; % Because values greater than 1 are due to rounding/the fact that N is calculated a posteriori rather than as part of the optimisation

p7c = plot(dataStore2S(1:8,1),NSvals(1:8),'Color',IBMpink,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p15c = plot(dataStore2S(9:16,1),NSvals(9:16),'Color',IBMorange,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p30c = plot(dataStore2S(17:24,1),NSvals(17:24),'Color',IBMyellow,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);
p60c = plot(dataStore2S(25:32,1),NSvals(25:32),'Color',IBMblue,'LineWidth',3,'LineStyle',':','Marker','x','MarkerSize',10);

xlabel('Sampling effort (E)');
ylabel('N_S^* (proportion of maximum possible value)');

leg = legend([p7c p15c p30c p60c],'\Delta = 7 days','\Delta = 15 days','\Delta = 30 days','\Delta = 60 days','NumColumns',2,'Orientation','vertical','Location','southwest');

ylim([0.75,1.05]);
% sampleSizes = (dataStore2S(1:32,1).*dataStore2S(1:32,2))./((dataStore2S(1:32,6)-dataStore2S(1:32,5)+1)/365)


%% -----------------------------------------------------------------------------------------------
% EXTRACT DATA

% Extract sample size and sample interval vectors
startDayVec = unique(dataStore2S(:,3))';
endDayVec = unique(dataStore2S(:,4))';

startSize = length(startDayVec);
endSize = length(endDayVec);

% Extract matrix of the baseline EDP
baselineEDP = reshape(dataStore2S(:,5),[endSize startSize]);
baselineEDP_smooth = imgaussfilt(baselineEDP,1);

% Extract matrix of the new EDP
newEDP = reshape(dataStore2S(:,6),[endSize startSize]);
newEDP_smooth = imgaussfilt(newEDP,1);

percChange = 100*(newEDP_smooth-baselineEDP_smooth)./baselineEDP_smooth;

optPoint = min(min(percChange));
[endMin,startMin]=find(percChange==optPoint);
optStart = startDayVec(startMin);
optEnd = endDayVec(endMin);