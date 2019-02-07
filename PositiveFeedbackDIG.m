% The beginning of a lumped parameter model of positive feedback
% Lumps transcription and translation together
% Separates protein A and B
% Assumes inducer is not used up in the reaction
% Adapted specifically for DIG system
% Output per cell
% Change to minutes

% Jeremy Primus
% March 1, 2017

% specify rate constants
A = 6.022e23;       % Avagadro's #
%kfi = 0.000005;     % induction rate constant
PA = 850e-13;       % permeability constant - 50*nm/s = 3um/min*area of plasmodesmata~283um converted to liters
kb = 0.00001;       % binding parameter - includes value for constituitive expression of unstable DIG-binder protein  
% kdi = 5e-9;       % actual DIG dissociation constant
kId = 2.5e-2;       % inducer degradation rate (half-life approx. 2 days)
kBd = 1.5e-3;       % B degradation rate constant (8 hrs) (half-life of luciferase ~ 2.5 hrs)

% hill equation parameters
n = 2;          % hill coefficient 
Vm = 100;    % vmax (max transcription rate)
kp = 1e-5;    % half max value

% specify species quantities - amount (could possibly do concentration)
Pa = 1;    % quantity of promoter A

% define the system of equations
f = @(t,x) [-kId*x(1);                                                    % dIo/dt - mass kinetics - outside [Ind]: degradation
             PA*(x(1) - x(2)) - kId*x(2);                                 % !!!dIi/dt -fick's law + mass kinetics - inside [Ind]: degradation
             kb*x(2) - kId*x(3);                                          % dA/dt - mass kinetics - # of DIG/protein complex 
             ((0.5*Vm*(x(3))^(n))/(kp + (x(3))^(n))) - kBd*x(4); 
             ((Vm*(x(4)+x(5))^n)/(kp + (x(4)+x(5))^n)) - kBd*x(5);];                   % dB/dt - induction by protein A + Hill equation

% specify initial condition
Io = 1e-6;    % Inducer concentration in mol/L
Ii = 0;       % Inducer concentration in the cell.  Probably always start at 0
Ao = 1;       % Inducer*uptake*binding.  In the absence of I, this can represent leaky expression terms for A
Bo = 1;       % and B
Ion = Io*A;   % Inducer concentration in molecules/L

% call on ode23 to solve coupled differential equations
[t,xx] = ode23(f,[0 5000],[Ion Ii Ao Bo 0]);

Luc = xx(:,4)+xx(:,5);
% plot the system dynamics
figure()
plot(t,xx(:,3:4))
title('Positive feedback')
xlabel('time (s)'), ylabel('Substance')
legend('Inducer','A','B','Location','best')
figure()
plot(t,xx(:,1))
xlabel('time (min)'), ylabel('Inducer oustide')
figure()
plot(t,xx(:,2))
xlabel('time (min)'), ylabel('Inducer inside')
figure
plot(t,xx(:,3))
xlabel('time (min)'), ylabel('Protein A')
figure
semilogy(t,xx(:,4))
xlabel('time (min)'), ylabel('Protein B')
figure
plot(t,xx(:,4))
xlabel('time (min)'), ylabel('Protein B')
plot(t,Luc(:),t,xx(:,4:5))
xlabel('time (min)'), ylabel('Protein')
legend('Luc Output','Due to A*','Due to B','Location','best')
