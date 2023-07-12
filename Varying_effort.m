%% Code for Fig 3C: computing the baseline EDP for a range of sample sizes and intervals
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code uses the functions 'runSpread_1' and 'runSampling_1' to compute the
% baseline EDP across a specified range of sample sizes and intervals. 
clear; close all;

%% ------------------------------------------------------------------------
% DEFINE MODEL PARAMETERS

% Population size (P=P_C)
P = 1000; 
% Transmission coefficient for 'Detectable' crops (beta_C)
r=0.05; beta = r/P;
% Transmission scaling factor for 'Undetectable' crops (epsilon_C)
epsilon = 0.015;
% Duration of crop 'Undetectable' period (gamma_C)
gamma = 452;
% Initial numbers of 'Detectable' and 'Undetectable' crops
D0 = 0; U0 = 1; 

numYears = 20;

% Vector of sample sizes to consider
effortVec = 0.5:0.5:4.5;
% Vector of sample intervals to consider
DminVec = 7:7:63;

% Number of spread simulations to perform
numSims = 1000; 
% Number of sampling simulations to perform
numRuns = numSims;
% Maximum run time for spread and sampling simulations
tFinal = 5000; 
% Specify whether to display progress messages
progress = "no";

%% ------------------------------------------------------------------------
% SPECIFY DISTRIBUTION OF START DAYS (COMMENT OUT AS APPROPRIATE)

% Seasonal transmission rate (according to fitted vector numbers)
C1 = 0.007265; a = 0.0004574; g0 = 0.064913; tem = 80;
day_vals = 1:1:365;
Pvec = C1*exp(g0*(day_vals-tem)-(a/2)*(day_vals-tem).^2);
omega = Pvec/mean(Pvec);
pdf = omega/sum(omega);
cdf = cumtrapz(pdf);
r = rand(1,numSims);

startDays = floor(interp1(cdf,day_vals,r)+1);
omega1 = omega;

% Seasonal symptom development rate (according to temperature)
temp = [9 9 11 14 19 23 26 26 22 18 13 10];
temp(temp<17)=0;
temp = temp-min(temp); temp = temp/mean(temp);
temp = interp1(1:1:12,temp,linspace(1,12,365));

omega2 = temp;

%% ------------------------------------------------------------------------
% END USER INPUT

%% ------------------------------------------------------------------------
% GENERATE EPIDEMIC CURVES

simData = runTimeSpread_VSDR_1(P,D0,U0,beta,epsilon,gamma,omega1,omega2,numYears,numSims,startDays,tFinal,progress);

%% ------------------------------------------------------------------------
% PERFORM SAMPLING SCHEME FOR RANGE OF SAMPLE SIZES AND INTERVALS

% Create empty matrix to store generated values
EDPstore = zeros(length(effortVec),length(DminVec));
EDTstore = zeros(length(effortVec),length(DminVec));
EDRstore = zeros(length(effortVec),length(DminVec));

for i=1:length(effortVec) % Iterate over sample sizes
    effort = effortVec(i)
    
    for j=1:length(DminVec) % Iterate over sample intervals
        Dmin = DminVec(j)
        
        % Print progress message if specified
        if progress == "yes"; fprintf(strcat('Sample size =',32,num2str(effort),',',32,'Sample interval =',32,num2str(Dmin),'\n')); end
        
        sampleSize = min(P,round(effort*Dmin));

        % Perform sampling simulations and store results
        [sampleData, EDP, EDT, EDR, perc95] = runWindowSampling_1(simData,numRuns,P,0,364,sampleSize,Dmin,tFinal,progress);
        EDPstore(i,j) = EDP;     
        EDTstore(i,j) = EDT; 
        EDRstore(i,j) = EDR; 
    end
end

%% ------------------------------------------------------------------------
% PLOT

figure();
levels = [0 1 2:2:20 20:5:55];
% Create contour plot of EDPs expressed as percentage of total population
[C,h] = contourf(100*EDPstore'/P,levels);
% clim([0 55]);
% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Baseline EDP (% of population)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Set x axis label and ticks
xlabel('Sampling effort (E)')
xticks = 1:1:length(effortVec);
xticklabels = strsplit(num2str(effortVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

% Set y axis label and ticks
ylabel('Sample interval (D_{min})')
yticks = 1:1:length(DminVec);
yticklabels = strsplit(num2str(DminVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

box off; set(gca,'Fontsize',16,'Linewidth',2);
