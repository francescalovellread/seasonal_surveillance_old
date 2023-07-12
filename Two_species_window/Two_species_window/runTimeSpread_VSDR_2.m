function simData = runTimeSpread_VSDR_2(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,omega1,omega2,numYears,numSims,startDays,tFinal,progress)

    % INPUT
    % popSize: total population size.
    % numSentinels: number of sentinels in the population
    % initialI: total initial number of symptomatic individuals (crops and sentinels)
    % initialC: total inital number of cryptic/asymptomatic individuals (crops and sentinels)
    % betaC: transmission coefficient for symptomatic crops
    % betaS: transmission coefficient for symptomatic sentinels
    % epsilonC: transmission coefficient scaling factor for cryptic/asymptomatic crops
    % epsilonS: transmission coefficient scaling factor for cryptic/asymptomatic sentinels
    % gammaC: cryptic/asymptomatic period for crops
    % gammaS: cryptic/asypmtomatic period for sentinels
    % omega1: 1x365 vector containing values of time dependent transmission factor for each day of the year
    % omega2: 1x365 vector containing values of time dependent symptom developement factor for each day of the year
    % numYears: number of years to initialise the system with
    % numSims: number of simulations to run
    % startDays: 1xnumSims  vector containing values of the start day for infection to arrive
    % tFinal: maximum time for each simulation to run
    % progress: specifies whether progress messages are displayed ("yes" or "no")
    
    % OUTPUT
    % sim_data: a (no_sims x 7) cell array containing simulation results. Each row corresponds
    % to a run; columns contain vectors for t, Ic (symptomatic crops), Is (symptomatic sentinels),
    % Cc (cryptic crops), Cs (cryptic sentinels), Sc (susceptible crops) and Ss (susceptible
    % sentinels).
    
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
    
    P = popSize; Ps = numSentinels; Pc = P-Ps; sentinel_prop = Ps/P;
    bc = betaC; bs = betaS; ec = epsilonC; es = epsilonS; gc = gammaC; gs = gammaS;
    
    % Create cell array for storing results
    simData = cell(numSims,7);
    
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
        startDay = startDays(i);
       
        % On each run, divide initial infected individuals between crops and sentinels
        % according to their relative proportions in the population
        I0 = initialI; C0 = initialC;
        Ic0 = 0; Is0 = 0; Cc0 = 0; Cs0 = 0; randsI = rand(1,I0); randsC = rand(1,C0);
        for j = 1:I0
            if randsI(j)<=sentinel_prop
                Is0 = Is0+1;
            else
                Ic0 = Ic0+1;
            end
        end
        for j = 1:C0
            if randsC(j)<=sentinel_prop
                Cs0 = Cs0+1;
            else
                Cc0 = Cc0+1;
            end
        end
        % Compute initial numbers of susceptible crops and sentinels
        Sc0 = Pc-Ic0-Cc0; Ss0 = Ps-Is0-Cs0;
    
        % Set initial conditions
        t = startDay; Ic = Ic0; Is = Is0; Cc = Cc0; Cs = Cs0; Sc = Sc0; Ss = Ss0;
        % Crop infection: 1. Sentinel infection: 2.
        % Crop symptoms: 3. Sentinel symptoms: 4.
        T1 = 0; T2 = 0; T3 = 0; T4 = 0;
        Q1 = log(1/rand(1)); Q2 = log(1/rand(1)); Q3 = log(1/rand(1)); Q4 = log(1/rand(1));
    
        % Create empty vectors to store results
        vecLength = 1+2*(Sc0+Ss0)+Cc0+Cs0;
%         tvec = zeros(1,vecLength);
%         Icvec = zeros(1,vecLength); Isvec = zeros(1,vecLength);
%         Ccvec = zeros(1,vecLength); Csvec = zeros(1,vecLength);
%         Scvec = zeros(1,vecLength); Ssvec = zeros(1,vecLength);
%         T1vec = zeros(1,vecLength); T2vec = zeros(1,vecLength); T3vec = zeros(1,vecLength); T4vec = zeros(1,vecLength); 
%         Q1vec = zeros(1,vecLength); Q2vec = zeros(1,vecLength); Q3vec = zeros(1,vecLength); Q4vec = zeros(1,vecLength); 
    
        mystore = zeros(15,vecLength);
        
        % Populate first entries with initial conditions
        index = 1; 
% 
        mycol = [t, Ic0, Is0, Cc0, Cs0, Sc0, Ss0, T1, T2, T3, T4, Q1, Q2, Q3, Q4];
        mystore(:,index) = mycol;
%         tvec(index)=t; Icvec(index)=Ic0; Isvec(index)=Is0; Ccvec(index)=Cc0; Csvec(index)=Cs0; Scvec(index)=Sc0; Ssvec(index)=Ss0;
%         T1vec(index) = T1; T2vec(index) = T2; T3vec(index) = T3; T4vec(index) = T4;
%         Q1vec(index) = Q1; Q2vec(index) = Q2; Q3vec(index) = Q3; Q4vec(index) = Q4;
    
    
        while (t<tFinal && Sc+Ss+Cc+Cs>0)
    
            % Compute individual reaction propensities
            a1 = Sc*(bc*Ic+bs*Is+ec*bc*Cc+es*bs*Cs);
            a2 = Ss*(bc*Ic+bs*Is+ec*bc*Cc+es*bs*Cs);
            a3 = (1/gc)*Cc;
            a4 = (1/gs)*Cs;
    
            % Compute constant 'right hand sides' of integral equation
            RHS1 = (Q1-T1)/a1;
            RHS2 = (Q2-T2)/a2;
            RHS3 = (Q3-T3)/a3;
            RHS4 = (Q4-T4)/a4;
    
            % Compute taus by solving integral equations
            tVal12 = interp1qr(days,ints1,t);

            ttauVal1 = tVal12+RHS1;
            ttau1 = interp1qr(ints1,days,ttauVal1);
            tau1 = ttau1 - t; if isnan(tau1); tau1 = Inf; end
    
            ttauVal2 = tVal12+RHS2;
            ttau2 = interp1qr(ints1,days,ttauVal2);
            tau2 = ttau2 - t; if isnan(tau2); tau2 = Inf; end
    
            tVal34 = interp1qr(days,ints2,t);

            ttauVal3 = tVal34+RHS3;
            ttau3 = interp1qr(ints2,days,ttauVal3);
            tau3 = ttau3 - t; if isnan(tau3); tau3 = Inf; end
    
            ttauVal4 = tVal34+RHS4;
            ttau4 = interp1qr(ints2,days,ttauVal4);
            tau4 = ttau4 - t; if isnan(tau4); tau4 = Inf; end
    
            % Find minimum tau to determine which reaction happens and update system state accordingly
            taus = [tau1 tau2 tau3 tau4];
            [~,I] = min(taus);
    
            if I==1 % then reaction 1 occurs (crop infection)
                tau = tau1;
                factor34 = interp1qr(days,ints2,t+tau)-tVal34;
                T1 = Q1;
                T2 = T2 + a2*RHS1;
%                 T3 = T3 + a3*(interp1qr(days_shift,ints2,t+tau)-tVal34);
%                 T4 = T4 + a4*(interp1qr(days_shift,ints2,t+tau)-tVal34);
                T3 = T3 + a3*factor34;
                T4 = T4 + a4*factor34;

                Q1 = Q1+log(1/rand(1));
                Sc = Sc-1; Cc = Cc+1;
            
            elseif I==2 % then reaction 2 occurs (sentinel infection)
                tau = tau2;
                factor34 = interp1qr(days,ints2,t+tau)-tVal34;
                T1 = T1 + a1*RHS2;
                T2 = Q2;
%                 T3 = T3 + a3*(interp1qr(days_shift,ints2,t+tau)-tVal34);
%                 T4 = T4 + a4*(interp1qr(days_shift,ints2,t+tau)-tVal34);
                T3 = T3 + a3*factor34;
                T4 = T4 + a4*factor34;

                Q2 = Q2+log(1/rand(1));
                Ss = Ss-1; Cs = Cs+1;

            elseif I==3 % then reaction 3 occurs (crop symptoms)
                tau = tau3;
                factor12 = interp1qr(days,ints1,t+tau)-tVal12;
%                 T1 = T1 + a1*(interp1qr(days_shift,ints1,t+tau)-tVal12);
%                 T2 = T2 + a2*(interp1qr(days_shift,ints1,t+tau)-tVal12);
                T1 = T1 + a1*factor12;
                T2 = T2 + a2*factor12;
                T3 = Q3;
                T4 = T4 + a4*RHS3;
                
                Q3 = Q3+log(1/rand(1));
                Cc = Cc-1; Ic = Ic+1;

            else % reaction 4 occurs (sentinel symptoms)
                tau = tau4;
                factor12 = interp1qr(days,ints1,t+tau)-tVal12;
%                 T1 = T1 + a1*(interp1qr(days_shift,ints1,t+tau)-tVal12);
%                 T2 = T2 + a2*(interp1qr(days_shift,ints1,t+tau)-tVal12);
                T1 = T1 + a1*factor12;
                T2 = T2 + a2*factor12;
                T3 = T3 + a3*RHS4;
                T4 = Q4;
                
                Q4 = Q4+log(1/rand(1));
                Cs = Cs-1; Is = Is+1;
            end
    
            % Update t and index and store new values
            t = t+tau;
            index = index + 1; 
            mycol = [t, Ic, Is, Cc, Cs, Sc, Ss, T1, T2, T3, T4, Q1, Q2, Q3, Q4];
            mystore(:,index) = mycol;

        end
    
        % Store results
        simData{i,1} = mystore(1,:);
        simData{i,2} = mystore(2,:); simData{i,3} = mystore(3,:);
        simData{i,4} = mystore(4,:); simData{i,5} = mystore(5,:);
        simData{i,6} = mystore(6,:); simData{i,7} = mystore(7,:);
    end
    
    elapsedTime = toc;
    if progress == "yes"
        fprintf(strcat('DONE! (',num2str(elapsedTime),32,'secs)\n',num2str(numSims),32,'incidence curves generated.\n\n'));
    end
end

