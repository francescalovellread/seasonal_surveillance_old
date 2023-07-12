%% Code for Fig 6: Optimising the total number of sentinels in the population and the sample
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code iterates over specified values of the sample size and sample interval,
% and in each instance applies a Bayesian optimisation algorithm to compute the
% optimal number of sentinels to include in the population and in the sample. It
% also computes the corresponding reduction in EDP compared to the baseline level.

function sentinels_RUN(ID)

%% ------------------------------------------------------------------------
% DEFINE MODEL PARAMETERS

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

% Vector of sample sizes to consider
sampleSizeVec = 25:25:200;
% Vector of sample intervals to consider
sampleIntervalVec = 30:15:150;

% Number of spread simulations to perform
numSims = 2500; 
% Number of sampling simulations to perform
numRuns = numSims;
% Maximum run time for spread and sampling simulations
tFinal = 8000; 
% Specify whether to display progress messages
progress = "no";

numYears = 20;

% Number of iterations for Bayesian optimisation algorithm
numIterations = 30;
% Specify number of parallel workers
nWorkers = 1;
% Define file path for save location
savePath = './results/';
% Define random number generator
rng('shuffle');


%% ------------------------------------------------------------------------
% SPECIFY DISTRIBUTION OF START DAYS according to FITTED vector numbers

C1 = 0.007265;
a = 0.0004574;
g0 = 0.064913;
tem = 80;

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

% omega1 = ones(1,365);
% omega2 = ones(1,365);
% startDays = ones(365);


%% ------------------------------------------------------------------------
% END USER INPUT

%% ------------------------------------------------------------------------
% EXTRACT SAMPLE SIZE AND SAMPLE INTERVAL FROM LISTS

[p,q] = meshgrid(sampleSizeVec,sampleIntervalVec);
pairs = [p(:) q(:)];
pairs_choose = pairs(ID,:);
sampleSize = pairs_choose(1);
sampleInterval = pairs_choose(2);

%% ------------------------------------------------------------------------
% COMPUTE BASELINE EDP

simData = runTimeSpread_VSDR_1(Pc,D0,U0,betaC,epsilonC,gammaC,omega1,omega2,numYears,numSims,startDays,tFinal,progress);
[~, EDP, ~, ~, ~] = runWindowSampling_1(simData,numRuns,Pc,0,364,sampleSize,sampleInterval,tFinal,progress);
baselineEDP = EDP;

%% ------------------------------------------------------------------------
% PERFORM BAYESIAN OPTIMISATION

% Define the number of sentinels in the population as an optimisable variable
totalNumSentinels = optimizableVariable('Stot',[0,350],'Type','integer');
% Define the number of sentinels in the sample as an optimisable variable        
numSentinelsSampled = optimizableVariable('Ssamp',[0,sampleSize],'Type','integer');

% Specify objective function input (function defined at end of script)
fun = @(z)objFunc(Pc+z.Stot,z.Stot,D0,U0,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,omega1,omega2,numYears,numSims,startDays,tFinal,progress,numRuns,sampleSize-z.Ssamp,z.Ssamp,sampleInterval,baselineEDP);
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

% Compute resultant EDP at optimum
resultantEDP = 0.01*valueAtOpt*baselineEDP+baselineEDP;

% Store results
results = [sampleSize, sampleInterval, optSentinelsAdded, optSentinelsSampled, valueAtOpt, baselineEDP, resultantEDP];

%% ------------------------------------------------------------------------
% WRITE RESULTS TO .TXT FILES

% Make new directory to store results
mkdir(savePath)

T = table(Pc,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,D0,U0,numSims,numRuns,tFinal,numIterations);
if ID==1
    writetable(T,[savePath 'params.txt'],'Delimiter','tab');
end

filename = sprintf('results_%d.txt',ID);
writematrix(results,[savePath filename],'Delimiter','tab');

end

