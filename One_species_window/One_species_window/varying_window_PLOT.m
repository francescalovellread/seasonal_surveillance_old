%% Code for Supplementary Figs S1-14: Sensitivity analyses
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% % This code takes in the output from 'sensitivity_RUN.m' and uses it to plot
% Supplementary Figs S1-14

%% -----------------------------------------------------------------------------------------------
% READ IN DATA

% Provide path to data folder
cd '/Users/francescalovell-read/Documents/DPhil/Chapter_4/CHAPTER_4_CODE/One_species_window/ARC/minD=15/effort=2/results';

% Make directory of all results files
FileList  = dir(fullfile(pwd,'results_*.txt'));

% Import and store data from results files
dataStore = [];
for ID = 1:numel(FileList)
    filename = sprintf('results_%d.txt',ID);
    data = importdata(filename);
    dataStore = [dataStore;data]; 
end

% Read in parameter table
params = readtable('params.txt');
P = params.P;
effort = params.effort;
minD = params.minD;

% -----------------------------------------------------------------------------------------------
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


% optStore = [effort, minD, optPoint, optStart, optEnd];
% 
% filename = sprintf('optResults_%.1f_%d.txt',effort,minD);
% 
% savePath = '/Users/francescalovell-read/Documents/DPhil/Chapter_4/CHAPTER_4_CODE/One_species_window/ARC/OPTIMAL_RESULTS/';
% 
% writematrix(optStore,[savePath filename],'Delimiter','tab');


%% ------------------------------------------------------------------------
% PLOT PERCENTAGE CHANGE

figure();

colormap(1-winter);
[C,h] = contourf(percChange);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, '% Reduction in EDP compared to USB', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

xlabel('Day sampling begins (D_{start})')
xticks = 1:2:length(startDayVec);
xticklabels = strsplit(num2str(startDayVec(1:2:end)));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

% Set y axis label and ticks
ylabel('Day sampling ends (D_{end})')
yticks = 1:2:length(endDayVec);
yticklabels = strsplit(num2str(endDayVec(1:2:end)));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

set(gca,'fontsize',16,'Linewidth',2); box off;

hold on;
mymin = plot(startMin,endMin);
mymin.Marker = 'square';
mymin.MarkerSize = 10;   
mymin.MarkerFaceColor = 'k';

% title('effort=0.5');

