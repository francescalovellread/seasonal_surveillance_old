

sampleSizeVec = 25:5:200;
sampleIntervalVec = 30:5:150;
[c,d]=meshgrid(sampleSizeVec,sampleIntervalVec);
effort = c./d;

levels = [0 0.5 1 1.5 2 2.5 3 3.5 4];
[C,h] = contourf(effort,levels);

clim([0 6]);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end


% % Define colourbar
colbar = colorbar;
ylabel(colbar, 'Sampling effort E = N/\Delta', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:5:length(sampleSizeVec);
xticklabels = strsplit(num2str(sampleSizeVec(xticks)));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:3:length(sampleIntervalVec);
yticklabels = strsplit(num2str(sampleIntervalVec(yticks)));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;