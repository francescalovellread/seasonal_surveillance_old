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

% Vector of start days to consider
startDvec = 0:10:160;
% Vector of end days to consider
endDvec = 200:10:360;

% Number of spread simulations to perform
numSims = 200; 
% Number of sampling simulations to perform
numRuns = numSims;
% Maximum run time for spread and sampling simulations
tFinal = 8000; 
% Specify whether to display progress messages
progress = "no";

numYears = 20;

sampleInterval = minD;


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
% GENERATE EPIDEMIC CURVES
simData = runTimeSpread_VSDR_1(P,D0,U0,beta,epsilon,gamma,omega1,omega2,numYears,numSims,startDays,tFinal,progress);

%% ------------------------------------------------------------------------
% RUN SAMPLING UNIFORMLY THROUGHOUT THE YEAR
startD = 0; endD = 364; prop = 1; 
sampleSize = min(P,round(effort*sampleInterval/prop));

[~, EDP, ~, ~, ~] = runWindowSampling_1(simData,numRuns,P,startD,endD,sampleSize,sampleInterval,tFinal,progress);
baselineEDP = EDP;

%% ------------------------------------------------------------------------
% RUN SAMPLING RESTRICTED TO SPECIFIED WINDOW
startD = pairs_choose(1); endD = pairs_choose(2); prop = (endD-(startD-1))/365;
sampleSize = min(P,round(effort*sampleInterval/prop));

[~, EDP, ~, ~, ~] = runWindowSampling_1(simData,numRuns,P,startD,endD,sampleSize,sampleInterval,tFinal,progress);
newEDP = EDP;

%% ------------------------------------------------------------------------
% STORE RESULTS
results = [effort, minD, startD, endD, baselineEDP, newEDP];

%% ------------------------------------------------------------------------
% WRITE RESULTS TO .TXT FILES

% Make new directory to store results
mkdir(savePath)

T = table(effort,minD,P,beta,epsilon,gamma,D0,U0,numSims,numRuns,tFinal,numYears,minD,effort);
if ID==1
    writetable(T,[savePath 'params.txt'],'Delimiter','tab');
end

filename = sprintf('results_%d.txt',ID);
writematrix(results,[savePath filename],'Delimiter','tab');

end





