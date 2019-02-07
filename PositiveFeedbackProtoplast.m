%% Positive Feedback DIG
% Recoding in preparation for ONR yearly report
% This code is intended to simulate the expected results of future
% protoplast experiment

% Jeremy Primus
% June 8, 2017

%% Parameters

A = 6.022e23;       % Avagadro's number
v = 1.41e-11;       % volume of a typical protoplast

% specify rate constants
% time units in minutes
d1 = 0.00012;   % inducer degradation rate (calculated from 4 day half-life)
d2 = 0;         % protein A degradation 
d3 = 0.0038;    % protein A* degradation (estimated as 2.5 hr half-life)- to be fit
d6 = 0.0038;    % protein B degradation (estimated as 2.5 hr half-life)- to be fit
d5 = 0.0038;    % luciferase degradation (estimated as 2.5 hr half-life)

kd = 6e-9;      % A* dissociation constant (Baker)




% hill equation parameters
n1 = 2;                 % Hill coefficient (best guess is 2:  dimerization of Gal4 to pUAS - Xie et al. 2000)
Vm = 30;              % vmax (max transcription rate) - mRNAs/min (TR for HIS3 gene)
kp = (50e-20);       % half max value (promoter half-saturation for transcription factor) (50 - 1000 nM - Yeast Genetic Networks): same for both Hill fns - same promoter

n2 = 2;
a = 0.6;

%% The system of equations
% define the system of equations
f = @(t,x) [-d1*x(1);                                                                                       % dI/dt - DIG - mass kinetics: - degradation
             0;                                                                                             % dA/dt - unbound A - mass kinetics: constituitive expression - degradation - binding + unbinding  
            -d3*x(3);                                                                              % dA*/dt - active bound A* - mass kinetics: binding - unbinding - degradation
             ((a*Vm*(x(3))^(n1))/(kp + (x(3))^(n1))) - d5*x(4);                                             % dC/dt - production of luciferase due to A* protein             
             ((Vm*(x(3))^(n1))/(kp + (x(3))^(n1))) + ((Vm*(x(6))^(n2)))/(kp + (x(6))^(n2)) - d5*x(5);    % dC/dt - total production of luciferase
             ((Vm*(x(6))^(n2))/(kp + (x(6))^(n2))) - d6*x(6)];                                            % dB/dt - positive feedback                      
         
%% Initial Conditions
% specify initial condition
Io = 0.1e-6;                         % Inducer concentration in mol/L
ao = 0.1e-9;                         % Steady state of protein A prior to induction
Ao = (Io*ao)/kd;                     % 
Bo = 0.1e-9;                         % and B
Co = 0.1e-9;                         % Leaky expression of C
Ic = [Io ao Ao Co Co Bo];            % initial conditions in concentration
IcN = (A/v)*[Io ao Ao Co Co Bo];     % initial conditions in number
%% Integration/Simulation
% call on ode23 to solve coupled differential equations
[t,xx] = ode45(f,[0 5000],Ic);

% plot the system dynamics
figure()
plot(t,xx(:,1:6))
title('Positive feedback')
xlabel('time (s)'), ylabel('Substance')
legend('Inducer','A','A*','C due to A*', 'Total C', 'Positive Feedback','Location','best')