%% Code for Supplementary Figs S1-14: Sensitivity analyses
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% % This code takes in the output from 'sensitivity_RUN.m' and uses it to plot
% Supplementary Figs S1-14

%% -----------------------------------------------------------------------------------------------
% READ IN DATA

% Provide path to data folder
cd '/Users/francescalovell-read/Documents/DPhil/Chapter_4/CHAPTER_4_CODE/Fig_6/full_seasonality/results';

% Make directory of all results files
FileList  = dir(fullfile(pwd,'results_*.txt'));

% %%
dataStore = [];
for ID = 1:540
    filename = sprintf('results_%d.txt',ID);
    if isfile(filename)
         % File exists.
    else
         fprintf(strcat(num2str(ID),'\n'));
    end
end
%%
% Import and store data from results files
dataStore = [];
for ID = 1:numel(FileList)
    filename = sprintf('results_%d.txt',ID);
    data = importdata(filename);
    dataStore = [dataStore;data]; 
end

% Read in parameter table
params = readtable('params.txt');
Pc = params.Pc;

%% -----------------------------------------------------------------------------------------------
% EXTRACT DATA

% Extract sample size and sample interval vectors
sampleSizeVec = unique(dataStore(:,1))';
sampleIntervalVec = unique(dataStore(:,2))';

nSize = length(sampleSizeVec);
nInterval = length(sampleIntervalVec);

% Extract matrix of the optimal total number of sentinels in the population
optTotSentinels = reshape(dataStore(:,3),[nInterval nSize]);

% Extract matrix of the optimal proportion of sentinels in the sample
optSentinelsSampled = reshape(dataStore(:,4),[nInterval nSize]);
optSentProp = optSentinelsSampled./min(optTotSentinels,sampleSizeVec);

% Extract matrix of the reduction in EDP at the optimum
EDPreduction = reshape(dataStore(:,5),[nInterval nSize]);

% Extract matrix of the resultant EDP at the optimum
EDPresultant = reshape(dataStore(:,7),[nInterval nSize]);
EDPresultant = 100*EDPresultant/Pc;

% Compute matrix of the relative sentinel utility
utility = -EDPreduction./optTotSentinels;

%% ------------------------------------------------------------------------
% PLOT OPTIMAL NUMBER OF SENTINELS

figure();

% Create contour plot of optimal number of sentinels
optTotSentinels_smooth = imgaussfilt(optTotSentinels,1);
levels = 0:20:350;
[C,h] = contourf(optTotSentinels_smooth,levels);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Optimal number of sentinels (P_S^*)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:5:nSize;
xticklabels = strsplit(num2str(sampleSizeVec(xticks)));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:3:nInterval;
yticklabels = strsplit(num2str(sampleIntervalVec(yticks)));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;

% hold on
% % pgon = polyshape([1.25 10 10 1.25],[6 6 15 15]);
% line1 = line([1.25 10],[6 6]);
% line2 = line([10 10],[6 15]);
% line3 = line([10 1.15],[15 15]);
% line4 = line([1.25 1.25],[15 6]);
% 
% line1.Color = 'r'; line1.LineWidth = 3;
% line2.Color = 'r'; line2.LineWidth = 3;
% line3.Color = 'r'; line3.LineWidth = 3;
% line4.Color = 'r'; line4.LineWidth = 3;
% 
% xlim([1,nSize]);
% ylim([1,13]);

%% ------------------------------------------------------------------------
% PLOT OPTIMAL SENTINEL PROPORTION

figure(); 

% Create contour plot of optimal sentinel proportion
levels = [0:0.2:0.8 0.9 0.95 0.99 1];
optSentProp_smooth = imgaussfilt(optSentProp,1);
[C,h] = contourf(optSentProp_smooth,levels);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = 'r'; end
clim([0 1]);

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Proportion of sentinels sampled at optimum', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:5:nSize;
xticklabels = strsplit(num2str(sampleSizeVec(xticks)));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:3:nInterval;
yticklabels = strsplit(num2str(sampleIntervalVec(yticks)));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;

%% ------------------------------------------------------------------------
% PLOT REDUCTION IN EDP

figure(); 

% Create contour plot of EDP reduction
EDPreduction_smooth = imgaussfilt(EDPreduction,1);
[C,h] = contourf(EDPreduction_smooth);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Change in EDP from baseline at optimum (%)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);
% colbar.TickLabels{end}='0^+';

% Format x axis
xlabel('Sample size (N)');
xticks = 1:5:nSize;
xticklabels = strsplit(num2str(sampleSizeVec(xticks)));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:3:nInterval;
yticklabels = strsplit(num2str(sampleIntervalVec(yticks)));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;

%% ------------------------------------------------------------------------
% PLOT RESULTANT EDP

figure();

% Create contour plot of resultant EDP
EDPresultant_smooth = imgaussfilt(EDPresultant,1);
levels = [0:1:10 10:2:20];
[C,h] = contourf(EDPresultant_smooth,levels);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'EDP at optimum (% of crop population)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:5:nSize;
xticklabels = strsplit(num2str(sampleSizeVec(xticks)));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:3:nInterval;
yticklabels = strsplit(num2str(sampleIntervalVec(yticks)));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;

%% ------------------------------------------------------------------------
% PLOT SENTINEL UTILITY

figure(); 

% Create contour plot of sentinel utility
utility_smooth = imgaussfilt(utility,1);
levels = [-1 0:0.1:5];
[C,h] = contourf(utility_smooth,levels);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end
clim([0 1.2]);

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Relative utility of sentinel', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:5:nSize;
xticklabels = strsplit(num2str(sampleSizeVec(xticks)));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:3:nInterval;
yticklabels = strsplit(num2str(sampleIntervalVec(yticks)));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;
