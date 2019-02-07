%% Kevin TATA Variants
%   Jeremy Primus
%   June 13, 2017

%% Description
% This general piece of code can be used for any type of dual luciferase
% assay where one type of luciferase measure is used for normalization

%% Initialization 
% test group list in the order it appears in the data set
Variants = {'Wildtype PNOS','PNOS Consensus','PNOS A5C','PNOS A4G', 'PNOS TC89GG','PNOS 1-4 complement','PNOS 5-8 complement','PNOS 4bp deletion'};
rep = 3;      % enter the number of replicates here 
tol = 0.06e9; % renilla tolerance.  Renilla RLU must be > tol to be included
stnd = 1;     % index of the reference standard in test group

%% Normalization and exclusion
Firefly = vec2mat(FLuc, rep);   % converts to matrix. Rows of test groups.  Columns of replicates
Renilla = vec2mat(RLuc, rep);   % converts to matrix. Rows of test groups.  Columns of replicates
Firefly(Renilla < tol) = NaN;   % excludes values without adaquate renilla for normalization

% normalize firefly output to constituitive renilla output
Fnorm = Firefly./Renilla;

% average the included replicates
Flucavg = sum(Fnorm, 2, 'omitnan')./(rep - sum(isnan(Fnorm),2));
% standard deviation of replicates 
dev = std(Fnorm,0,2,'omitnan')/;

% scale the firefly matrix
ref = Flucavg(stnd);
RelPromoters = Flucavg/ref;
dev2 = dev/ref;

% sort data by strength - Strongest to weakest
[B,K] = sort(RelPromoters);
K = K';
K2 = fliplr(K);
B = B';
B2 = fliplr(B);

% create a bar graph of relative promoter strengths
figure
b = bar(B2);
c = b.FaceColor;
b.FaceColor = [0.9 0.3 0.5];
title('PNOS Variants')
set(gca,'xticklabel',Variants(K2));
set(gca,'xticklabelrotation',90);
hold on
errorbar(1:8, B2, dev2(K2),'.')


