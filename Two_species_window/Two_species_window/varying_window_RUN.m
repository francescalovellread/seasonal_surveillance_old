%% Code for Fig 3C: computing the baseline EDP for a range of sample sizes and intervals
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code uses the functions 'runSpread_1' and 'runSampling_1' to compute the
% baseline EDP across a specified range of sample sizes and intervals. 

function varying_window_RUN(ID)

%% ------------------------------------------------------------------------
% DEFINE MODEL PARAMETERS


effort = 1;
minD = 30; 


% Crop population size (PC)
Pc = 1000; 
% Transmission coefficient for 'Detectable' crops (betaC) 
r=0.05; betaC = r/Pc;
% Transmission coefficient for 'Detectable' sentinels (betaS)
betaS = betaC;
% Transmission scaling factor for 'Undetectable' crops (epsilonC)
epsilonC = 0.015;
% Transmission scaling factor for 'Undetectable' sentinels (epsilonS)
epsilonS = 0.1;
% Duration of crop 'Undetectable' period (gammaC)
gammaC = 452;
% Duration of sentinel 'Undetectable' period (gammaS)
gammaS = 49;
% Initial numbers of 'Detectable' and 'Undetectable' plants
D0 = 0; U0 = 1;

% Vector of start days to consider
startDvec = 0:10:160;
% Vector of end days to consider
endDvec = 200:10:360;

% Number of spread simulations to perform
numSims = 1000; 
% Number of sampling simulations to perform
numRuns = numSims;
% Maximum run time for spread and sampling simulations
tFinal = 8000; 
% Specify whether to display progress messages
progress = "no";

numYears = 20;
sampleInterval = minD;

% Number of iterations for Bayesian optimisation algorithm
numIterations = 30;
% Specify number of parallel workers
nWorkers = 1;
% Define file path for save location
savePath = './results/';
% Define random number generator
rng('shuffle');

%%
% Seasonal transmission rate (according to fitted vector numbers)
C1 = 0.007265; a = 0.0004574; g0 = 0.064913; tem = 80;

day_vals = 1:1:365;
Pvec = C1*exp(g0*(day_vals-tem)-(a/2)*(day_vals-tem).^2);
omega = Pvec/mean(Pvec);
pdf = omega/sum(omega);
cdf = cumtrapz(pdf);
r = rand(1,numSims);

startDays = floor(interp1(cdf,day_vals,r)+1);

% Seasonal symptom development rate (according to temperature)
temp = [9 9 11 14 19 23 26 26 22 18 13 10];
temp(temp<17)=0;
temp = temp-min(temp); temp = temp/mean(temp);
temp = interp1(1:1:12,temp,linspace(1,12,365));

omega1 = omega;
omega2 = temp;

%% ------------------------------------------------------------------------
% END USER INPUT

%% ------------------------------------------------------------------------
% EXTRACT START DAY AND END DAY FROM LISTS

[p,q] = meshgrid(startDvec,endDvec);
pairs = [p(:) q(:)];
pairs_choose = pairs(ID,:);

%% ------------------------------------------------------------------------
% BASELINE CASE WITH UNIFORM SAMPLING

simData = runTimeSpread_VSDR_1(Pc,D0,U0,betaC,epsilonC,gammaC,omega1,omega2,numYears,numSims,startDays,tFinal,progress);
startD = 0; endD = 364; prop = 1; 
sampleSize = min(Pc,round(effort*sampleInterval/prop));

[~, EDP, ~, ~, ~] = runWindowSampling_1(simData,numRuns,Pc,startD,endD,sampleSize,sampleInterval,tFinal,progress);
baselineEDP = EDP;

%% ------------------------------------------------------------------------
% COMPUTE PROPORTION OF YEAR AND SAMPLE SIZE
startD = pairs_choose(1);
endD = pairs_choose(2);
prop = (endD-startD)/365;
sampleSize = min(Pc,round(effort*sampleInterval/prop));

%% ------------------------------------------------------------------------
% PERFORM BAYESIAN OPTIMISATION

% Define the number of sentinels in the population as an optimisable variable
totalNumSentinels = optimizableVariable('Stot',[0,350],'Type','integer');
% Define the number of sentinels in the sample as an optimisable variable        
numSentinelsSampled = optimizableVariable('Ssamp',[0,sampleSize],'Type','integer');

% Specify objective function input (function defined at end of script)
fun = @(z)objFunc(Pc+z.Stot,z.Stot,D0,U0,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,omega1,omega2,numYears,numSims,startDays,tFinal,progress,numRuns,sampleSize-z.Ssamp,z.Ssamp,sampleInterval,startD,endD);
% Perform Bayesian optimisation
optResults = bayesopt(fun,[totalNumSentinels,numSentinelsSampled],'IsObjectiveDeterministic',false,'MaxObjectiveEvaluations',numIterations,'AcquisitionFunctionName','expected-improvement-plus','XConstraintFcn',@objConstraint,'UseParallel',false,'PlotFcn',[]);

%% ------------------------------------------------------------------------
% EXTRACT AND STORE RESULTS

% Extract optimal parameters
optimalParams = optResults.XAtMinEstimatedObjective;
% Extract optimal number of sentinels to include in population
optSentinelsAdded = optimalParams.Stot;
% Extract optimal number of sentinels to include in sample
optSentinelsSampled = optimalParams.Ssamp;

% Extract value of objective function at optimum
valueAtOpt = optResults.MinEstimatedObjective;

% Store results
results = [effort, minD, startD, endD, baselineEDP, optSentinelsAdded, optSentinelsSampled, valueAtOpt];

%% ------------------------------------------------------------------------
% WRITE RESULTS TO .TXT FILES

% Make new directory to store results
mkdir(savePath)

T = table(effort,minD,Pc,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,D0,U0,numSims,numRuns,tFinal,numIterations,numYears);
if ID==1
    writetable(T,[savePath 'params.txt'],'Delimiter','tab');
end

filename = sprintf('results_%d.txt',ID);
writematrix(results,[savePath filename],'Delimiter','tab');

end

% %
% %
% figure()
% % contourf(store2')
% % levels = [0 30:1:50 50:5:100];
% 
% % levels = [-50:2:0 0:5:50];
% adv = 100*(storeee-bl)/bl;
% adv = imgaussfilt(adv,1)';
% [C,h] = contourf(adv);
% % 
% % levels = [0 30:1:50 50:5:100];
% % contourf(imgaussfilt(storeee,1)',levels);
% 
% % % Set contour and label styles and specify manual labelling
% h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
% mylabels = clabel(C,h,'manual','color','k','FontSize',16);
% for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end
% 
% colormap(1-winter);
% % Set x axis label and ticks
% xlabel('Sampling begins')
% xticks = 1:1:length(startDvec);
% xticklabels = strsplit(num2str(startDvec));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
% 
% % Set y axis label and ticks
% ylabel('Sampling ends')
% yticks = 1:1:length(endDvec);
% yticklabels = strsplit(num2str(endDvec));
% set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
% 
% box off; set(gca,'Fontsize',16,'Linewidth',2);
% colorbar
% 
% 
% minimum = min(min(adv));
% [y,x]=find(adv==minimum);
% hold on
% mymin = plot(x,y);
% mymin.Marker = 'square';
% mymin.MarkerSize = 10;   
% mymin.MarkerFaceColor = 'k';






