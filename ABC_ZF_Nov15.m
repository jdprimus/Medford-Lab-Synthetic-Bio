%% Approximate Bayesian Computation (ABC)
% November 28, 2018
% Jeremy Primus

%% Updates
% Update November 28th: 
% Run which began on November 27th appears to be complete (see 'ABC_ZF_Nov15.mat').
% Analyze results, put notes below:
% 'A' values may not appear to be saving correctly.  Suspicious of the
% number of 1s in the results.
% Made some changes to the order in which A's were saved,testing run:
% 'ABC_ZF_Nov27.mat'
% This may be alright, as the 1s only show up in the later distributions.
% Updating comments and checking code for obvious mistakes, as this will be the last run of this software in the foreseeable future. 

%% Description
% Approximate Bayesian Computation (ABC), in an implemenation of a
% Sequential Markov Chain (SMC), similar to as described in DOI:
% 10.1002/env.2353 (provided by Jennifer Hoeting) is used to estimate the
% parameters in a mechanistic model to describe a plant-based genetic
% circuit.

% Dataset:
% Written for plate induction experiment conducted on (date?).  Raw data at
% (location) was processed with 'Shell3.m' to produce quantified data at
% C:\Users\Medford Lab\Desktop\QuantifiedZFData. Processing notes can be
% found here as well.  The quantified data was compiled using
% 'Post_Analysis.m'

% authored dependancies:
% compileExperimentData.m
% zfMechanisticModelNonNeg
% rpert.m

% Parent: ABC_ZF_meanasdistancefcn
% Modified code to implement the changes discussed in September 7th meeting
% with Ashok Prasad.
% Promising results were obtained from parent: see ABCmeanasdistancefcnAug30.mat  
% Changes implemented from parent:
% uniform priors
% save 'A' values for data model

% Notes:
% This code has not been thoroughly tested in its entirety, and the
% correctness and effect of the updates have not been tested at all.
%
% These things should be checked:
% assess summary statistics
%      Take a more thorough look at data distributions (which appear to be lognormal)
%      and ensure that they justify the choice of summary statistics
% assess distance function
%      Given the choice of summary statistics, is the chosen distance
%      function (one measure for the entire dataset) the most appropriate?
% tune step size via 'shape' parameter for optimal acceptance - (test via summit)
%      More for performance than influence on results.  A perturbation
%      kernel would also be a more elegant way to step.
% is there a better way to sample 'A'? 
%       The 'quick and dirty' implementation used here should work, and at
%       this point does not justify a change unless the parameter distribtuions
%       suggest one.
% validate and justify the tolerance cascade
%       The tolerance was selected based upon trial and error.  However,
%       the tolerance should be set by some sort of variance in the true
%       data set.
% do distributions need to be proper?
%       Nothing has been done to validate whether the distributions used
%       and obtained are proper, although this may not be necessary.
% particle weighting
%       The change to uniform distributions should have obviated some of
%       the reasons for this change.  However, the best form of particle
%       weighting should correct the final distibutions to account for how
%       frequently parameters were drawn, and how close the simulation lies
%       to the 'truth' (via summary statistics and distance function).
% sample paths, acceptance rate, and GR statistic
%       Of course, with any Bayesian method these should be assessed and
%       reported.

% This code could be applied to a several datasets, such as the liquid and
% soil induction experiments conducted in September 2018 (partially
% quantified).  The following recommendations could be implemented for
% added modularity to different datasets:
    % extract plants_per_treat, timepoints, treatments from experimental data
    % ensure dataset comes in order of timepoints (ex. T1-T5)
    % put in explicit commands to remove timepoints(ex. T0)
    % separate parameter model, so this code does not depend on the exact form
    % of the mechanistic model
    % extract tolerance from the variability in the true data
    % automate the process of tuning step size
    % generalize code
    % optimize code

%% Parameter initialization
clear all;
% path = '/home/jdprimus@colostate.edu/QuantifiedZFData';   % The path to experimental data - this is the summit path
path = 'u:\QuantifiedZFData';                               % The path to experimental data
data = compileExperimentData(path);                         % compiles the experimental data to fit
LogMeans = readtable('ZFconexpLogMeansOrdered.xlsx');
data = [data LogMeans];
plants_per_treat = max(data.PlantNo_);                      % plants per treatment
tMapMin = 60*[17.176 24 41.17 48 65];                       % timepoints of data collection
timepoints = length(tMapMin);                               % no. of timepoints
nA = 6.022e23;                                              % Avagadro's number
treatmentName = {'Control' '0.1' '1' '5' '20' '50'};
treatment = [0 0.1 1 5 20 50]*1e-6*nA;                      % inducer treatment groups

% ABC parameters
tol_fb = [3.2 2.9 2.5 2.2 2.0 1.8]; 
N = 10000;                                   % initial number of steps per t    
T = plants_per_treat;                        % number of intermediates
shape = 4000;                                % tuning parameter for sampling from rpert

% ODE solver options
interrupt_time = 30;                                                 % sets max solve time for a parameter set in the ODE solver
outputFun= @(t,y,flag)interruptFun(t,y,flag,interrupt_time);
opts = odeset('NonNegative', 1:5, 'OutputFcn', outputFun); 

% initialize output parameters
posterior_theta = {};                       % matrix to save posterior distributions
posterior_A_da = {};                        % matrix to save direct activation scalar
posterior_A_fb = {};                        % matrix to save positive feedback scalar
theta_star = ones(1,15,T);                  % matrix to save all proposed parameters
t = 1;
% parpool(24) - summit
parpool(48)

%% Run ABC
while t <= T 
    accepted = false(N,1);                      % matrix of logical values to indicate whether the proposed parameters meet acceptance criteria
    parfor n = 1:N       % parameter draws iterable
        if t == 1
            %% Prior distributions - change to uniform priors?
            % time unit - minutes
            % 90% CI given after each prior

            % deltaI - degradation rate constant of inducer - prior based on a 28.5
            % minute half life for atmospheric oxidation of OHT suggested by: http://www.chemspider.com/Chemical-Structure.4447687.html
            deltaI_star = random('Uniform', 0.00204963, 0.112928); % 0.00204963-0.112928 min^-1 = 6-338 min half-life

            % kf - forward rate constant of inducer/NEV binding.  From information at: https://www.annualreviews.org/doi/full/10.1146/annurev-biophys-070816-033639?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub%3Dpubmed&
            kf_star = 60*10^(random('Uniform', 3.4, 6.6))/nA;     % 60*10^3.4 - 60*10^6.6 (M*min)^-1 = 4.2e-21 - 6.6e-18 (molecules*min)^-1

            % kb - backward rate constant of inducer/NEV dissociation from KD given at: https://www.sciencedirect.com/science/article/pii/S0303720798001300?via%3Dihub
            kb_star = nA*kf_star*0.16e-9;

            % K - constituitive NEV production of 35S promoter
            K_star = random('Uniform', 2.5, 149.8); % 2.5 - 149.8 molecules/min

            % deltaN - rate constant of NEV degradation
            deltaN_star = random('Uniform', 0.00204963, 0.112928); % 0.00204963-0.112928 min^-1 = 6-338 min half-life

            % alphaZ - leaky expression of zinc finger production
            alphaZ_star = random('Uniform', 0.1, 6); % 0.1 - 6 molecules/min

            % Vm1 - maximal production rate of zinc finger from 10XN1 promoter.  Used
            % information on maximal transcription and translation rates from: http://book.bionumbers.org/what-is-faster-transcription-or-translation/
            Vm1_star = random('Uniform', 2.5, 149.8); % 2.5 - 149.8 molecules/min

            % n1 - Hill coefficient for the 10XN1 promoter
            n1_star = random('Uniform', 1.6, 12.5); % 1.6 - 12.5

            % Kp1 - half max saturation concentration for the 10XN1 promoter
            Kp1_star = 10^random('Uniform', 0, 30);

            % deltaZ - rate constant for degradation of the zinc finger
            deltaZ_star = random('Uniform', 0.00204963, 0.112928); % 0.00204963-0.112928 min^-1 = 6-338 min half-life

            % alphaL - leaky expression of luciferase
            alphaL_star = random('Uniform', 0.1, 6); % 0.1 - 6 molecules/min

            % Vm2 - maximal production rate of luciferase/zinc finger from the pUAS promoter
            Vm2_star = random('Uniform', 2.5, 149.8); % 2.5 - 149.8 molecules/min

            % n2 - Hill coefficient for the pUAS promoter
            n2_star = random('Uniform', 1.6, 12.5); % 1.6 - 12.5

            % Kp2 - half max saturation concentration for the pUAS promoter
            Kp2_star = 10^random('Uniform', 0, 30);

            % deltaL - degradation rate constant of luciferase
            deltaL_star = random('Uniform', 0.00204963, 0.112928); % 0.00204963-0.112928 min^-1 = 6-338 min half-life
        else
            %% Intermediate Proposal distributions - change to kernel perturbation
            % time unit - minutes
            % 90% CI given after each proposal
            intermediate_theta = cell2mat(posterior_theta(t-1));
            nonzeroindex = intermediate_theta(:,1)~=0;
            intermediate_theta = intermediate_theta(nonzeroindex,:);
            h = size(intermediate_theta,1);
            
            % draw an index for each parameter from the posterior
            draw = random('Discrete Uniform', h);

            % deltaI - degradation rate constant of inducer - prior based on a 28.5
            % minute half life for atmospheric oxidation of OHT suggested by: http://www.chemspider.com/Chemical-Structure.4447687.html
            deltaI_star = rpert(0, 1, min(intermediate_theta(draw,1)), shape); % 0.00125106-0.0730666 min^-1 = 9.5-554 min half-life

            % kf - forward rate constant of inducer/NEV binding.  From information at: https://www.annualreviews.org/doi/full/10.1146/annurev-biophys-070816-033639?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub%3Dpubmed&
            kf_o = log10(intermediate_theta(draw,2)*nA/60);   % invert draw from posterior_params to a useable hyperparameter
            kf_star = 60*10^(random('Normal', kf_o, 1))/nA;     % 60*10^3.4 - 60*10^6.6 (M*min)^-1 = 4.2e-21 - 6.6e-18 (molecules*min)^-1

            % kb - backward rate constant of inducer/NEV dissociation from KD given at: https://www.sciencedirect.com/science/article/pii/S0303720798001300?via%3Dihub
            kb_star = nA*kf_star*0.16e-9;

            % K - constituitive NEV production of 35S promoter
            K_star = rpert(0, 200, min(intermediate_theta(draw,4),200), shape); % 2.5 - 149.8 molecules/min

            % deltaN - rate constant of NEV degradation
            deltaN_star = rpert(0, 1, min(intermediate_theta(draw,5),1), shape); % 0.00125106-0.0730666 min^-1 = 9.5-554 min half-life

            % alphaZ - leaky expression of zinc finger production
            alphaZ_star = rpert(0, 100, min(intermediate_theta(draw,6),100), shape); 

            % Vm1 - maximal production rate of zinc finger from 10XN1 promoter.  Used
            % information on maximal transcription and translation rates from: http://book.bionumbers.org/what-is-faster-transcription-or-translation/
            Vm1_star = rpert(0, 200, min(intermediate_theta(draw,7),200), shape); % 2.5 - 149.8 molecules/min

            % n1 - Hill coefficient for the 10XN1 promoter
            n1_star = rpert(0, 20, min(intermediate_theta(draw,8),20), shape); % 1.6 - 12.5

            % Kp1 - half max saturation concentration for the 10XN1 promoter
            Kp1_star = random('LogNormal', log(intermediate_theta(draw,9)), 4);

            % deltaZ - rate constant for degradation of the zinc finger
            deltaZ_star = rpert(0, 1, min(intermediate_theta(draw,10),1), shape); % 0.00125106-0.0730666 min^-1 = 9.5-554 min half-life

            % alphaL - leaky expression of luciferase
            alphaL_star = rpert(0, 100, min(intermediate_theta(draw,11),100), shape); 

            % Vm2 - maximal production rate of luciferase/zinc finger from the pUAS promoter
            Vm2_star = rpert(0, 200, min(intermediate_theta(draw,12),200), shape); % 2.5 - 149.8 molecules/min

            % n2 - Hill coefficient for the pUAS promoter
            n2_star = rpert(0, 20, min(intermediate_theta(draw,13),20), shape); % 1.6 - 12.5

            % Kp2 - half max saturation concentration for the pUAS promoter
            Kp2_star = random('LogNormal', log(intermediate_theta(draw,14)), 4);

            % deltaL - degradation rate constant of luciferase
            deltaL_star = rpert(0, 1, min(intermediate_theta(draw,15),1), shape); % 0.00125106-0.0730666 min^-1 = 9.5-554 min half-life
        end
        %% save all proposed parameters to track sample paths
        theta_star(n,:,t) = [deltaI_star kf_star kb_star K_star deltaN_star alphaZ_star...
                       Vm1_star n1_star Kp1_star deltaZ_star alphaL_star Vm2_star n2_star Kp2_star deltaL_star];
        %% Pass parameters to Process model
        % generate predicted initial conditions from parameters
        % implement parfeval when all plants and treatments are added -
        % does this help if all workers are already slotted for the parfor

        treatment_accepted = false(length(treatment),1);            % intialize boolean variable indicating whether the parameter set was accepted for all treatments
        noOfTreatments = length(treatment);
        A_da_star(:,:,n) = zeros(noOfTreatments, plants_per_treat);
        A_fb_star(:,:,n) = zeros(noOfTreatments, plants_per_treat);
        

        da_exp_mean = ones(noOfTreatments,5);
        fb_exp_mean = ones(noOfTreatments,5);
        da_exp_dev = ones(noOfTreatments,5);
        fb_exp_dev = ones(noOfTreatments,5);
        da_sim_mean = ones(noOfTreatments,5);
        fb_sim_mean = ones(noOfTreatments,5);
        da_sim_dev = ones(noOfTreatments,5);
        fb_sim_dev = ones(noOfTreatments,5);
        scale_temp_da = ones(noOfTreatments, plants_per_treat);
        scale_temp_fb = ones(noOfTreatments, plants_per_treat);

        
        for i = 1:noOfTreatments                                % loop over treatments
            treatmentTable = data(strcmp(data.Treatment, treatmentName{i}),:);                                % obtain all data for treatment

            IC = [treatment(i) K_star/deltaN_star 0 alphaZ_star/deltaZ_star alphaL_star/deltaL_star;...
                  treatment(i) K_star/deltaN_star 0 alphaZ_star/deltaZ_star alphaL_star/deltaL_star];         % intial conditions

            % simulate molecular outputs from process model
            [gx, gt, g_index, fx, ft, f_index] = zfMechanisticModelNonNeg(deltaI_star, kf_star, kb_star, K_star, deltaN_star, alphaZ_star, Vm1_star, n1_star, ...
                                                            Kp1_star, Vm2_star, n2_star, Kp2_star, deltaZ_star, alphaL_star, deltaL_star, IC, opts, tMapMin);
            % check if integration timed out
            if isnan(sum(sum(gx))) == 0               % skip to next loop iteration if gx is NaN 
                % extract simulated data at timepoints nearest data collection timepoints                                                
                y_da_star = log(real(gx(g_index,5)));  % direct activation simulated data
                y_fb_star = log(real(fx(f_index,5)));  % positive feedback simulated data 
                %% pass the process model output to the data model
                
                % initialize empty arrays for loop over plants
                da_plant_data = ones(timepoints+1, plants_per_treat);
                fb_plant_data = ones(timepoints+1, plants_per_treat);                
                for j = 1:plants_per_treat         % loop over all plants
                    % direct activation 
                    da_plant_data(:,j) = treatmentTable.LogMean(strcmp(treatmentTable.Line, 'MSA103-13') & treatmentTable.PlantNo_ == j);    % extract relevant data
                    da_scale = da_plant_data(2:6,j)./y_da_star;                   % assess a scale factor to translate model output to each data point
                    scale_temp_da(i,j) = real(mean(da_scale));                             % determine mean scale factor for each direct activation plant
                    
                    % positive feedback
                    fb_plant_data(:,j) = treatmentTable.LogMean(strcmp(treatmentTable.Line, 'MSA107-18') & treatmentTable.PlantNo_ == j);
                    fb_scale = fb_plant_data(2:6,j)./y_fb_star;                   % assess a scale factor to translate model output to each data point
                    scale_temp_fb(i,j) = real(mean(fb_scale));                             % determine mean scale factor for each direct activation plant

%                                         figure()
%                     plot(tMapMin, fb_plant_data(2:6, j), tMapMin, A_fb(i,j)*y_fb_star)
                end

                
                A_fb_temp = scale_temp_fb;
                A_da_temp = scale_temp_da;
                
                % for this treatment: (both positive fb and da)
                % calculate mean of the true data
                da_exp_mean(i,:) = mean(da_plant_data(2:6,:),2)';
                fb_exp_mean(i,:) = mean(fb_plant_data(2:6,:),2)';
                % calculate standard deviation of the true data 
                da_exp_dev(i,:) = std(da_plant_data(2:6,:),0,2)';
                fb_exp_dev(i,:) = std(fb_plant_data(2:6,:),0,2)';
                % calculate mean of the simulated data
                da_sim_mean(i,:) = mean([A_da_temp(i,1)*y_da_star A_da_temp(i,2)*y_da_star A_da_temp(i,3)*y_da_star A_da_temp(i,4)*y_da_star A_da_temp(i,5)*y_da_star A_da_temp(i,6)*y_da_star]');
                fb_sim_mean(i,:) = mean([A_fb_temp(i,1)*y_fb_star A_fb_temp(i,2)*y_fb_star A_fb_temp(i,3)*y_fb_star A_fb_temp(i,4)*y_fb_star A_fb_temp(i,5)*y_fb_star A_fb_temp(i,6)*y_fb_star]');
                % calculate the std dev of the simulated data
                da_sim_dev(i,:) = std([A_da_temp(i,1)*y_da_star A_da_temp(i,2)*y_da_star A_da_temp(i,3)*y_da_star A_da_temp(i,4)*y_da_star A_da_temp(i,5)*y_da_star A_da_temp(i,6)*y_da_star]');
                fb_sim_dev(i,:) = std([A_fb_temp(i,1)*y_fb_star A_fb_temp(i,2)*y_fb_star A_fb_temp(i,3)*y_fb_star A_fb_temp(i,4)*y_fb_star A_fb_temp(i,5)*y_fb_star A_fb_temp(i,6)*y_fb_star]');

            end               
        end
        A_da_star(:,:,n) = scale_temp_da;
        A_fb_star(:,:,n) = scale_temp_fb;
        
% distance function            
    err_da(n,:,:) = da_exp_mean-da_sim_mean;
    distancefcn_da(n) = sqrt(sum(sum(err_da(n,:,:).^2)));
    err_fb(n,:,:) = fb_exp_mean-fb_sim_mean;
    distancefcn_fb(n) = sqrt(sum(sum(err_fb(n,:,:).^2)));
    if distancefcn_fb(n) < tol_fb(t)
        accepted(n) = 1
    end
    end
        %% criteria for moving on to next intermediate
        if sum(accepted) > 0                                                               % ensures some parameters have been accepted before moving on
            save('ABC_ZF_Nov27.mat')
            if isempty(posterior_theta) || length(posterior_theta) < t
                posterior_theta{t} = theta_star(accepted == 1,:,t);                        % move first working parameter sets to posterior/intermediate
                posterior_A_da{t} = A_da_star(:,:,accepted == 1);
                posterior_A_fb{t} = A_fb_star(:,:,accepted == 1);
            else
                posterior_A_da{t} = cat(3,posterior_A_da{t}, A_da_star(:,:,accepted == 1));  % concatonate additional working sets to posterior/intermediate
                posterior_A_fb{t} = cat(3,posterior_A_fb{t}, A_fb_star(:,:,accepted == 1));  % concatonate additional working sets to posterior/intermediate
                posterior_theta{t} = [posterior_theta{t}; theta_star(accepted == 1,:,t)];
            end
            if sum(posterior_theta{t}(:,1)~=0,1) >= 1000                                    % refine tolerance if 1000 parameter sets have been accepted at the current tolerance
                t = t+1;
                save('ABC_ZF_Nov28.mat')
            end
        end
end