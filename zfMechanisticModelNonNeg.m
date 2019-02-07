%% ZF Mechanistic Model - ALWAYS NONNegative
% Jeremy Primus
% July 9, 2018

function [gx, gt, g_index, fx, ft, f_index] = zfMechanisticModelNonNeg(deltaI, kf, kb, K, deltaN, alphaZ, Vm1, n1, Kp1, Vm2, n2, Kp2, deltaZ, alphaL, deltaL, IC, opts, tMapMin)
% define the system of equations
% x(1) = Inducer
% x(2) = NEV
% x(3) = active NEV
% x(4) = zinc finger
% x(5) = firefly luciferase

% parameter definitions
% deltaI - degradation rate constant of inducer
% kf - forward rate constant of inducer/NEV binding
% kb - backward rate constant of inducer/NEV dissociation
% K - constituitive NEV production of 35S promoter
% deltaN - rate constant of NEV degradation
% alphaZ - leaky expression of zinc finger production
% Vm1 - maximal production rate of zinc finger from 10XN1 promoter
% n1 - Hill coefficient for the 10XN1 promoter
% Kp1 - half max saturation concentration for the 10XN1 promoter
% deltaZ - rate constant for degradation of the zinc finger
% alphaL - leaky expression of luciferase
% Vm2 - maximal production rate of luciferase/zinc finger from the pUAS promoter
% n2 - Hill coefficient for the pUAS promoter
% Kp2 - half max saturation concentration for the pUAS promoter
% deltaL - degradation rate constant of luciferase

% direct activation model
g = @(t,x) [-deltaI*x(1) - kf*x(2)*x(1) + kb*x(3);                                              % inducer concentration              
             K - deltaN*x(2) - kf*x(2)*x(1) + kb*x(3);                                          % NEV concentration
             -deltaN*x(3) + kf*x(2)*x(1) - kb*x(3);                                             % active NEV concentration
             alphaZ + (Vm1*x(3)^n1)/(Kp1 + x(3)^n1) - deltaZ*x(4);                              % zinc finger concentration
             alphaL + (Vm2*x(4)^n2)/(Kp2 + x(4)^n2) - deltaL*x(5);];                            % firefly luciferase concentration
             
                try
                     [gt,gx] = ode23t(g,[0 4000],IC(1,:), opts);
                catch ME
                    if strcmp(ME.identifier,'interruptFun:Interrupt')
                       disp(ME.message);
                       [gx, gt, g_index, fx, ft, f_index] = deal(NaN);
                       return
                    % Do other things
                    else
                    rethrow(ME); % It's possible the error was due to something else
                    end
                end
            [~, g_index] = min(abs(gt-tMapMin));            % find index of timepoints nearest data collection timepoints

% positive feedback model
f = @(t,x) [-deltaI*x(1) - kf*x(2)*x(1) + kb*x(3);                                                      % inducer concentration              
             K - deltaN*x(2) - kf*x(2)*x(1) + kb*x(3);                                                  % NEV concentration
             -deltaN*x(3) + kf*x(2)*x(1) - kb*x(3);                                                     % active NEV concentration
             alphaZ + (Vm1*x(3)^n1)/(Kp1 + x(3)^n1) + (Vm2*x(4)^n2)/(Kp2 + x(4)^n2) - deltaZ*x(4);      % zinc finger concentration
             alphaL + (Vm2*x(4)^n2)/(Kp2 + x(4)^n2) - deltaL*x(5);];                                    % firefly luciferase concentration
                
                try
                    [ft,fx] = ode23t(f,[0 4000],IC(2,:), opts); 
                catch ME
                    if strcmp(ME.identifier,'interruptFun:Interrupt')
                       disp(ME.message);
                       [gx, gt, g_index, fx, ft, f_index] = deal(NaN);
                       return
                    % Do other things
                    else
                    rethrow(ME); % It's possible the error was due to something else
                    end
                end
            [~, f_index] = min(abs(ft-tMapMin));            % find index of timepoints nearest data collection timepoints

% plot the system dynamics
% figure()
% plot(gt,gx(:,1))
% figure()
% plot(gt,gx(:,2))
% figure()
% plot(gt,gx(:,3))
% figure()
% plot(gt,gx(:,4))
% figure()
% plot(gt,gx(:,5))
% title('Direct Activation')
% xlabel('time (s)'), ylabel('Luciferase')
% figure()
% plot(ft,fx(:,5))
% title('Positive feedback')
% xlabel('time (s)'), ylabel('Luciferase')

end