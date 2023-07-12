%% ------------------------------------------------------------------------
% DEFINE OBJECTIVE FUNCTION TO BE MINIMISED BY BAYESOPT

function OUTPUT = objFunc(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,omega1,omega2,numYears,numSims,startDays,tFinal,progress,numRuns,cropSampleSize,sentinelSampleSize,sampleInterval,openWindow,closeWindow)
    if numSentinels==0
        simData = runTimeSpread_VSDR_1(popSize,initialI,initialC,betaC,epsilonC,gammaC,omega1,omega2,numYears,numSims,startDays,tFinal,progress);
        [~, EDP, ~, ~, ~] = runWindowSampling_1(simData,numRuns,popSize,openWindow,closeWindow,cropSampleSize,sampleInterval,tFinal,progress);
        OUTPUT = EDP;
    else
        simData = runTimeSpread_VSDR_2(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,omega1,omega2,numYears,numSims,startDays,tFinal,progress);
        [~, ~, ~, ~, ECDP, ~, ~] = runWindowSampling_2(simData,numRuns,popSize,numSentinels,openWindow,closeWindow,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,progress);
        OUTPUT = ECDP;
    end
end
