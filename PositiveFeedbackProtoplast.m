%% Positive Feedback DIG
% Jeremy Primus
% June 8, 2017

%%  Description 
% Recoding in preparation for ONR yearly report
% This code is intended to simulate the expected results of future
% protoplast experiment

%%  Notes
% June 16, 2017
% Model performs well for concentrations, not for number
% Intuitively stable states seem very high in concentration units
% Need to get better parameters -- primarily diffision limited reaction
% rate
%% Parameters
clear all
A = 6.022e23;       % Avagadro's number
v = 1.41e-11;       % average volume of a protoplast

% specify rate constants
% time units in minutes
d1 = 0.00012;   % inducer degradation rate (calculated from 4 day half-life)
d2 = 0.13;      % protein A degradation  (5 min half-life Matalpha2 tag)
d3 = 0.0026;    % protein A* degradation (estimated as 2.5 hr half-life)- to be fit
d4 = 0.0046;    % luciferase degradation (estimated as 2.5 hr half-life)
d5 = 0.0046;    % luciferase degradation (estimated as 2.5 hr half-life)

kd = 6e-9;      % A* dissociation constant (Baker)
kf = 1e7;%*A*v; % estimated as the diffusion-limited rate constant (Holde 2002)
kb = kd*kf;
Ka = 0.03;


% Hill equation parameters
n1 = 2;                 % Hill coefficient for induction (best guess is 2:  dimerization of Gal4 to pUAS - Xie et al. 2000)
Vm = 0.2;                 % vmax (max transcription rate) - mRNAs/min (TR for HIS3 gene)
kp = (500e-3);%*A*v;    % half max value (promoter half-saturation for transcription factor) (50 - 1000 nM - Yeast Genetic Networks): same for both Hill fns - same promoter

n2 = 1.6;                 % Hill coefficient for positive feedback
a = 0.6;                % scale factor to compare VP16 to VP64 activation domain

% Explicit leaky expression 
L = 0.1e-5;

%% Alternate Parameter Set
clear all
A = 6.022e23;       % Avagadro's number
v = 1.41e-11;       % average volume of a protoplast

% specify rate constants
% time units in minutes
d1 = 0.00012;   % inducer degradation rate (calculated from 4 day half-life)
d2 = 0.13;      % protein A degradation  (5 min half-life Matalpha2 tag)
d3 = 0.0026;    % protein A* degradation (estimated as 2.5 hr half-life)- to be fit
d4 = 0.0046;    % luciferase degradation (estimated as 2.5 hr half-life)
d5 = 0.026;    % luciferase degradation (estimated as 2.5 hr half-life)

kd = 6e-9;      % A* dissociation constant (Baker)
kf = 1e7;%*A*v; % estimated as the diffusion-limited rate constant (Holde 2002)
kb = kd*kf;
Ka = 0.03;


% Hill equation parameters
n1 = 2;                 % Hill coefficient for induction (best guess is 2:  dimerization of Gal4 to pUAS - Xie et al. 2000)
Vm = 0.2;                 % vmax (max transcription rate) - mRNAs/min (TR for HIS3 gene)
kp = (500e-12);%*A*v;    % half max value (promoter half-saturation for transcription factor) (50 - 1000 nM - Yeast Genetic Networks): same for both Hill fns - same promoter

n2 = 3.74;                 % Hill coefficient for positive feedback
a = 0.6;                % scale factor to compare VP16 to VP64 activation domain

% Explicit leaky expression 
L = 0.00001;
%% The system of equations
% define the system of equations
f = @(t,x) [ -kf*x(2)*x(1) + kb*x(3) - d1*x(1);                                                             % dI/dt - DIG - mass kinetics: - degradation
             Ka - d2*x(2) - kf*x(2)*x(1) + kb*x(3);                                                         % dA/dt - unbound A - mass kinetics: constituitive expression - degradation - binding + unbinding  
             - d3*x(3) + kf*x(2)*x(1) - kb*x(3);                                                            % dA*/dt - active bound A* - mass kinetics: binding - unbinding - degradation
             ((a*Vm*(x(3))^(n1))/(kp + (x(3))^(n1))) - d4*x(4);                                             % dC/dt - production of luciferase due to A* protein             
             L + ((a*Vm*(x(3))^(n1))/(kp + (x(3))^(n1))) + ((Vm*(x(5))^(n2)))/(kp + (x(5))^(n2)) - d5*x(5);];   % dC/dt - total production of luciferase
                      
%% Initial Conditions
% specify initial condition
Io = 0.1e-6;                         % Inducer concentration in mol/L
ao = Ka/d2;                           % Steady state of protein A prior to induction
Ao = (Io*ao)/kd;                      % steady state quickly established for A*
Bo = 0;                               % B leaky expression   
Co = 0;                               % Leaky expression of C
Ic = [Io ao Ao Co Co];                % initial conditions in concentration
IcN = A*v*[Io ao Ao Co Co];           % initial conditions in number

%% Integration/Simulation
% call on ode23 to solve coupled differential equations
opts = odeset('NonNegative', 1:5);
[t,xx] = ode15s(f,[0 5000],Ic, opts);     % simulation time in minutes

%% Plot the system dynamics
% Plot all
figure()
plot(t,xx(:,:))
hold on
pfb = xx(:,5) - xx(:,4);            % Positive feedback luciferase
plot(t,pfb);
title('Positive feedback')
xlabel('time (min)'), ylabel('Substance')
legend('I','A','A*','Luciferase due to A*', 'Total Luciferase', 'Positive Feedback Luciferase','Location','best')       % legend for all outputs

% Plot only luciferase outputs
figure()
plot(t,xx(:,4:5))
hold on
pfb = xx(:,5) - xx(:,4);            % Positive feedback luciferase
plot(t,pfb);
title('Positive feedback')
xlabel('time (min)'), ylabel('Substance')
legend('Luciferase due to A*', 'Total Luciferase', 'Positive Feedback Luciferase','Location','best')                    % legend for Luciferase only output

%% Plot only MLB101 output
%figure()
hold on
plot(t,xx(:,5))
title('Stable States')
xlabel('time (min)'), ylabel('Substance')
%legend('No Inducer','Location','best')                    % legend for Luciferase only output
