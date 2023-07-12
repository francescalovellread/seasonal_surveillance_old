function [sampleData, EDP, EDT, EDR, perc95] = runWindowSampling_1(simData,numRuns,popSize,openWindow,closeWindow,sampleSizeIn,sampleIntervalIn,tFinal,progress)

% INPUT
% simData: cell array of simulated incidence curves, generated from runSpread_1.m
% numRuns: number of sampling runs to perform (must be <= number of simulations available
% in simData!)
% popSize: total population size
% sampleSize: number of plants to sample on each sampling round
% sampleInterval: sampling interval (time between each sampling round)
% tFinal: final permissible sample time
% progress: specifies whether progress messages are displayed ("yes" or "no")
% theta_S: sensitivity of detection during summer
% theta_W: sensitivity of detection during winter

% OUTPUT
% sampleData: a three-column matrix in which each row corresponds to a single
% sampling run. Entries in the first column are (total) discovery prevalences; entries in the
% second column are the corresponding discovery times; entries in the third column are the
% numbers of sampling rounds that occurred.ED
% EDP: the simulated expected discovery prevalence
% EDT: the simulated expected discovery time
% EDR: the simulated expected number of sampling rounds
% perc95: the 95th percentile of the detection prevalence

rng('shuffle');

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

timerSample1 = tic;

P = popSize;
N = sampleSizeIn; 
% Nout = sampleSizeOut;
Delta = sampleIntervalIn; 
% DeltaOut = sampleIntervalOut; 

% Extract number of simulations from simData
dataSize = size(simData); numSims = dataSize(1);
% If numRuns>numSims, return error message and exit
if numRuns > numSims
    fprintf(strcat('Error: number of sampling runs exceeds available number of simulations. Please set numRuns<=',num2str(numSims),'.\n\n'));
    return; 
end
selectedRuns = randsample(numSims,numRuns);
selectedSimData = simData(selectedRuns,:);
% Create empty matrix to store sampling data
sampleData = zeros(numRuns,3);
% Generate initial sampling times for each run. Initial sampling times are
% uniformly distributed on the interval [0,Delta] to mimic pathogen
% introduction at a random time relative to the sampling scheme.

%%******
initSampleTimes = rand(numRuns,1)*Delta;
%%******

if progress == "yes"
    fprintf('Running sampling simulations...\t')
end

% Run sampling simulations
for i=1:numRuns
%     initSampleTimes = rand(numRuns,1)*10
    sampleTime = initSampleTimes(i);
    detectionIndicator = 0;
    samplingRound = 1;

    if (openWindow<=mod(sampleTime,365) && mod(sampleTime,365)<closeWindow) % if we're in the 'window' period
        wInd = 1; % window indicator is ON
    else 
        wInd = 0; % window indicator is OFF
    end

    while (detectionIndicator == 0 && sampleTime <= tFinal)

%         if wInd == 1
%             N = Nin; Delta = DeltaIn;
%         else
%             N = Nout; Delta = DeltaOut;
%         end
        if wInd == 1
%             fprintf('SAMPLE!');
            sampleIndex = sum(selectedSimData{i,1} <= sampleTime); % Determine which time point on the inidence curve corresponds to the sample time
            if sampleIndex == 0
                sampleTime = sampleTime + Delta; % Move on to next sample time
                samplingRound = samplingRound + 1;
            else
                numDetectable = selectedSimData{i,2}(sampleIndex); % Compute total number of 'Detectable' plants
                numUndetectable = selectedSimData{i,3}(sampleIndex); % Compute total number of 'Undetectable' plants
                
%                 if (90<=mod(sampleTime,365) && mod(sampleTime,365)<273) % if we're in the 'summer' period (April-September)
%                     theta = theta_S; % detection sensitivity takes its 'summer' value
%                 else % we're in the 'winter' period (October-March)
%                     theta = theta_W; % detection sensitivity takes its 'winter' value
%                 end

                theta = symptoms(sampleTime);
    
                stateVec = rand(1,numDetectable); % Create vector containing a U[0,1] number for every 'Detectable' plant
                stateVec(stateVec<1-theta)=0; % Replace every entry less than detection threshold with a 0
                stateVec(stateVec>=1-theta)=1; % Replace every entry greater than detection threshold with a 1
                stateVec(P) = 0; % Add in 0s to represent the rest of the population that isn't Detectable 
         
                selectVec = randsample(stateVec,N); % Take a random sample of size N (without replacement) from this vector
                detSample = sum(selectVec); % Total number of 'Detectable' plants in the sample
                if detSample > 0 % Then a 'Detectable' plant has been sampled
                    detectionIndicator = 1;
                else
                    sampleTime = sampleTime + Delta; % Move on to next sample time
                    samplingRound = samplingRound + 1;
                    
                end
            end

        else
%             fprintf('SKIP!');
            sampleTime = sampleTime + Delta;
           
        end
        

        if (openWindow<=mod(sampleTime,365) && mod(sampleTime,365)<closeWindow) % if we're in the 'window' period
            wIndNew = 1; % window indicator is ON
        else 
            wIndNew = 0; % window indicator is OFF
        end


%         if wIndNew-wInd == 1 % then we've moved from winter to summer since last time!
%             % in that case: change sample size and sample interval and choose new sampling
%             % time at beginning of that window
% %             wIndNew-wInd
%             lastSampleTime = sampleTime-Delta;
% 
%             while (mod(lastSampleTime,365)<=openWindow || mod(lastSampleTime,365)>=closeWindow)
%                 lastSampleTime = lastSampleTime + Delta;
%             end
%             sampleTime = lastSampleTime;
% 
%         elseif wIndNew-wInd == -1 % then we've moved from summer to winter since last time!
%             % in that case: change sample size and sample interval and choose new sampling
%             % time at beginning of that window
% %             wIndNew-wInd
%             lastSampleTime = sampleTime-Delta;
% 
%             while (openWindow<=mod(lastSampleTime,365) && mod(lastSampleTime,365)<=closeWindow)
%                 lastSampleTime = lastSampleTime + Delta;
%             end
%             sampleTime = lastSampleTime;
% 
%         % else, we are still in the same window and nothing needs to change
%         end

        wInd = wIndNew;
%         pause
%         sampleTime
%         mod(sampleTime,365)
        


    end

    if detectionIndicator == 0
        sampleData(i,1) = P;
    else
        sampleData(i,1) = numDetectable + numUndetectable;
    end
    sampleData(i,2) = sampleTime;
    sampleData(i,3) = samplingRound;
end

EDP = mean(sampleData(:,1)); % Expected total discovery prevalence ('Detectable' and 'Undetectable')
EDT = mean(sampleData(:,2)); % Expected discovery time
EDR = mean(sampleData(:,3)); % Expected number of sampling rounds
perc95 = prctile(sampleData(:,1),95);

elapsedTimeSample1 = toc(timerSample1);
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTimeSample1),32,'secs)\n\n'));
end
end


% function theta = symptoms(t)
% 
%     precip = [46.7 49 45.3 33.9 22.0 13.1 10.8 16.5 42.7 61.6 70.1 58.0];
%     precip=precip-min(precip);
%     precip=precip/max(precip);
%     dryness = 1-precip;
%     
%     month = ceil(mod(t,365)*12/365);
%     if month==0; month=12; end
%     
%     theta = dryness(month);
% 
% end

% function theta = symptoms(t)
%     
%     theta = 1;
% 
% end
