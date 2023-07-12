%% Code for Fig 3C: computing the baseline EDP for a range of sample sizes and intervals
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code uses the functions 'runSpread_1' and 'runSampling_1' to compute the
% baseline EDP across a specified range of sample sizes and intervals. 

function varying_effort_RUN(ID)

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

% Vector of efforts to consider
effortVec = 0.5:0.5:4;
% Vector of sample intervals to consider
DminVec = 7:7:63;

% Number of spread simulations to perform
numSims = 10000; 
% Number of sampling simulations to perform
numRuns = numSims;
% Maximum run time for spread and sampling simulations
tFinal = 8000; 
% Specify whether to display progress messages
progress = "no";
% Define file path for save location
savePath = './results/';
% Define random number generator
rng('shuffle');

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
% EXTRACT SAMPLE SIZE AND SAMPLE INTERVAL FROM LISTS

[p,q] = meshgrid(effortVec,DminVec);
pairs = [p(:) q(:)];
pairs_choose = pairs(ID,:);
effort = pairs_choose(1);
Dmin = pairs_choose(2);

%% ------------------------------------------------------------------------
% GENERATE EPIDEMIC CURVES

simData = runTimeSpread_VSDR_1(P,D0,U0,beta,epsilon,gamma,omega1,omega2,numYears,numSims,startDays,tFinal,progress);

%% ------------------------------------------------------------------------
% PERFORM SAMPLING SCHEME FOR RANGE OF SAMPLE SIZES AND INTERVALS

sampleSize = min(P,round(effort*Dmin));
[~, EDP, ~, ~, ~] = runWindowSampling_1(simData,numRuns,P,0,364,sampleSize,Dmin,tFinal,progress);
results = [effort, Dmin, EDP];
    
%% ------------------------------------------------------------------------
% WRITE RESULTS TO .TXT FILES

% Make new directory to store results
mkdir(savePath)

T = table(P,beta,epsilon,gamma,D0,U0,numSims,numRuns,tFinal);
if ID==1
    writetable(T,[savePath 'params.txt'],'Delimiter','tab');
end

filename = sprintf('results_%d.txt',ID);
writematrix(results,[savePath filename],'Delimiter','tab');

end