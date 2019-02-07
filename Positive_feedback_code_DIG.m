% The beginning of a lumped parameter model of positive feedback
% Separates transcription and translation
% considers A-mRNA and B-mRNA to lead to the same protein
% Jeremy Primus
% Feb 20, 2017

% specify rate constants 
%kfi = 0.000005;     % induction rate constant
h = 0.3;            % uptake parameter
kdi = 9.625e-6;    % inducer degradation rate
kmAd = 0.00115;    % A-mRNA degradation rate constant
kmBd = 0.00115;    % B-mRNA degradation rate constant
kGd = 9.625e-5;    % protein of interest degradation rate constant
kt = 0.0125;      % translation rate constant

% hill equation parameters
n = 2;        % hill coefficient 
Vm = 0.5;       % vmax (max transcription rate)
kp = 10000000;   % half max value

% specify species quantities
Pa = 1; % quantity of promoter A

% define the system of equations
f = @(t,x) [ %- kfi*Pa*x(1);                 % dI/dt - mass kinetics - assumes inducer is used up in the reaction
             -kdi*x(1);                      % dI/dt - mass kinetics -
    h*x(1)*Pa - kmAd*x(2);                  
    (Vm*x(4)^n)/(kp + x(4)^n) - kmBd*x(3);   % dmB/dt - Hill equation
    kt*x(2) + kt*x(3) - kGd*x(4)];           % dG/dt - mass kinetics

% specify initial conditions
Io = 10;  % inducer molecules supplied
mAo = 0;    % these can be considered leaky expression terms for A-mRNA
mBo = 0;    % and B-mRNA
Go = 0;     % protein of interest initial amount (leave at 0 for now)

% call on ode23 to solve coupled differential equations
[t,xx] = ode23(f,[0 1000000],[Io mAo mBo Go]);


% plot the system dynamics
figure()
plot(t,xx(:,:))
title('Positive feedback')
xlabel('time (s)'), ylabel('Substance')
legend('Inducer','A mRNA','B mRNA','Protein of Interest','Location','best')