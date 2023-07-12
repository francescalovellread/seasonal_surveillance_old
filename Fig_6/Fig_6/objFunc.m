%% ------------------------------------------------------------------------
% DEFINE OBJECTIVE FUNCTION TO BE MINIMISED BY BAYESOPT

function OUTPUT = objFunc(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,omega1,omega2,numYears,numSims,startDays,tFinal,progress,numRuns,cropSampleSize,sentinelSampleSize,sampleInterval,baselineEDP)
    if numSentinels==0
        OUTPUT = 0;
    else
        simData = runTimeSpread_VSDR_2(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,omega1,omega2,numYears,numSims,startDays,tFinal,progress);
        [~, ~, ~, ~, ECDP, ~, ~] = runWindowSampling_2(simData,numRuns,popSize,numSentinels,0,364,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,progress);
        OUTPUT = 100*(ECDP-baselineEDP)/baselineEDP;
    end
end
