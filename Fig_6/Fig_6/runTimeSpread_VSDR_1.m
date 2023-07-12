function simData = runTimeSpread_VSDR_1(popSize,initialI,initialC,beta,epsilon,gamma,omega1,omega2,numYears,numSims,startDays,tFinal,progress)

    % INPUT
    % popSize: total population size.
    % initialI: initial number of symptomatic infected individuals
    % initialC: initial number of cryptic/asymptomatic infected individuals
    % beta: transmission coefficient for symptomatic individuals
    % epsilon: transmission coefficient scaling factor for cryptic/asymptomatic individuals
    % gamma: duration of asymptomatic period
    % omega1: 1x365 vector containing values of time dependent transmission factor for each day of the year
    % omega2: 1x365 vector containing values of time dependent symptom developement factor for each day of the year
    % numYears: number of years to initialise the system with
    % numSims: number of simulations to run
    % startDays: 1xnumSims  vector containing values of the start day for infection to arrive
    % tFinal: maximum time for each simulation to run
    % progress: specifies whether progress messages are displayed ("yes" or "no")
    
    % OUTPUT
    % sim_data: a (no_sims x 4) cell array containing simulation results. Each
    % row corresponds to a run; first, second, third and fourth columns contain
    % vectors for t, I, C and S respectively.
    
    if (progress ~= "yes" && progress ~= "no")
        fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
    end
    if (length(omega1) ~= 365 || length(omega2) ~= 365)
        fprintf('ERROR: "omega" vectors must contain 365 values\n\n'); return
    end
    if (length(startDays) ~= numSims)
        fprintf('ERROR: Vector "startDays" must contain numSims values\n\n'); return
    end

    tic
    
    P = popSize; I0 = initialI; C0 = initialC; b = beta; e = epsilon; g = gamma;
    % Compute initial number of susceptible individuals
    S0 = P-I0-C0;

    % Create cell array for storing results
    simData = cell(numSims,4);
   
    % Define days over which seasonal transmission is computed
    days = 0:1:numYears*365-1; days = days';
    % Repeat omega1 periodically over required number of years
    omega1_rep = repmat(omega1,1,numYears);
    % Integrate omega1 cumulatively over required number of years
    omega1_int = cumtrapz(omega1_rep); ints1 = omega1_int';
    % Repeat omega2 periodically over required number of years
    omega2_rep = repmat(omega2,1,numYears);
    % Integrate omega2 cumulatively over required number of years
    omega2_int = cumtrapz(omega2_rep); ints2 = omega2_int';

    if progress == "yes"
        fprintf('Running spread simulations...\t')
    end

    % Run stochastic simulations
    for i=1:numSims
        if mod(i,100)==0;
            i
        end
        startDay = startDays(i);
        % Set intial conditions
        t = startDay; I = I0; C = C0; S = S0; 
        % Infections: 1. Symptoms: 2
        T1 = 0; T2 = 0;
        Q1 = log(1/rand(1)); Q2 = log(1/rand(1));
    
        % Create empty vectors to store results
        vecLength = 1+2*S0+C0;
        mystore = zeros(8,vecLength);
        
        % Populate first entries with initial conditions
        index = 1; 

        mycol = [t, I0, C0, S0, T1, T2, Q1, Q2];
        mystore(:,index) = mycol;

        while (t<tFinal && (S+C)>0)

            % Compute individual reaction propensities
            a1 = b*S*(I+e*C); a2 = (1/g)*C; 
           
            % Compute constant 'right hand sides' of integral equation
            RHS1 = (Q1-T1)/a1;
            RHS2 = (Q2-T2)/a2;
    
            % Compute taus by solving integral equations
            tVal1 = interp1qr(days,ints1,t);
            ttauVal1 = tVal1+RHS1;
            ttau1 = interp1qr(ints1,days,ttauVal1);
            tau1 = ttau1 - t; if isnan(tau1); tau1 = Inf; end
            
            tVal2 = interp1qr(days,ints2,t);
            ttauVal2 = tVal2+RHS2;
            ttau2 = interp1qr(ints2,days,ttauVal2);
            tau2 = ttau2 - t; if isnan(tau2); tau2 = Inf; end
            
            % Find minimum tau to determine which reaction happens and update system state accordingly
            taus = [tau1 tau2];
            [~,min_tau] = min(taus);
% 
            if min_tau==1 % then reaction 1 occurs (infection)
                tau = tau1;

%                 T1 = T1 + a1*(interp1qr(days,ints1,t+tau)-interp1qr(days,ints1,t));
%                 T1 = T1 + a1*(interp1qr(days,ints1,t+tau)-tVal1);
                T1 = Q1;

%                 T2 = T2 + a2*(interp1qr(days,ints2,t+tau)-interp1qr(days,ints2,t));
                T2 = T2 + a2*(interp1qr(days,ints2,t+tau)-tVal2);

                Q1 = Q1+log(1/rand(1));
                C = C+1; S = S-1;

            else % reaction 2 occurs (symptom development)
                tau = tau2;

%                 T1 = T1 + a1*(interp1(days,ints1,t+tau)-interp1(days,ints1,t));
                T1 = T1 + a1*(interp1qr(days,ints1,t+tau)-tVal1);

%                 T2 = T2 + a2*(interp1(days,ints2,t+tau)-interp1(days,ints2,t));
%                 T2 = T2 + a2*(interp1qr(days,ints2,t+tau)-tVal2);
                T2 = Q2;
   
                Q2 = Q2+log(1/rand(1));
                I = I+1; C = C-1;

            end
           
            
            % Update t and index and store new values
            t = t+tau;
            index = index + 1; 
            mycol = [t, I, C, S, T1, T2, Q1, Q2];
            mystore(:,index) = mycol;

        end

        % Store results
        simData{i,1} = mystore(1,:);
        simData{i,2} = mystore(2,:);
        simData{i,3} = mystore(3,:);
        simData{i,4} = mystore(4,:);
    end
    
    elapsedTime = toc;
    if progress == "yes"
        fprintf(strcat('DONE! (',num2str(elapsedTime),32,'secs)\n',num2str(numSims),32,'incidence curves generated.\n\n'));
    end
end
