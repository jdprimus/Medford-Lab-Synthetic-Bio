%% Positive Feedback Histogram
% Jeremy Primus
% June 27, 2017

%% Description 
% Script to automate import  of quantified luciferase data and compile
% into pixel intensity histogram

%% Initialization
clear all;
cpath = 'C:\Users\Medford Lab\Downloads\OneDrive_1_6-20-2017\MLB101-SC_Control\Result Summary\HistogramData';  % path to control dataset
ipath = 'C:\Users\Medford Lab\Downloads\OneDrive_1_6-20-2017\MLB101-SC_Induced\Result Summary\HistogramData';  % path to induced dataset

tpi_portion = 'T001*_Whole*';  % designates the timepoint and portion of plant to histogram for the initial histogram
tpf_portion = 'T014*_Whole*';  % designates the timepoint and portion of plant to histogram for the final histogram

con_exc_plant = 19;             % plant to exclude from control group - set NaN to exclude none
ind_exc_plant = 8;              % plant to exclude from control group - set NaN to exclude none

bins = 100;                     % number of bins for histograms
logbins = logspace(1.5,6,bins); % sets histogram range
loggedcounts = zeros(1,bins);   % initializes bin counts to 0

line = 13;                      % set line to plot histograms for MLB101-12,13, or both (1213)

%% Control Initial
% histograms initial control distribution

files = dir(fullfile(cpath, tpi_portion));  % extract relevant files
if line == 12
    a = 1;
    b = length(files)/2;
elseif line == 13
    a = length(files)/2 + 1;
    b = length(files);
else
    a = 1;
    b = length(files);
end

% compile files into single table
for i=a:b
    filepixelmap = readtable(fullfile(cpath, files(i).name));
    [m,n] = size(filepixelmap);
    filepixelmap.Plant = (i+10)*ones(m,1);                                % append plant number as column to table
    if i == a
        totalpixelmap = filepixelmap;
    else    
        totalpixelmap = [totalpixelmap; filepixelmap];
    end
end

% table into arrays
pm = table2array(totalpixelmap);
keys = pm(:,1);         % assigns pixel intensities
values = pm(:,2);       % assigns counts at each intensity

% excludes a plant from histogram data
values(totalpixelmap.Plant == con_exc_plant ) = 0;

% compile a list of all unique pixel intensities
ukey = unique(keys);

% sum counts of pixel intensities 
for n = 1:length(ukey)
    totalcounts(n) = sum(values(keys == ukey(n)));
end

% bin pixel intensity counts logarithmically 
for j = 1:(length(logbins)-1)
    loggedcounts(j) = sum(totalcounts((ukey>logbins(j)) & (ukey<logbins(j+1))));
end
% finds the peak of the distribution
[~,I] = max(loggedcounts);
cipeak = logbins(I);

% create histogram
figure()
hci = histogram('BinEdges',logbins,'BinCounts', loggedcounts(1:bins-1),'Normalization','Probability','EdgeAlpha',0.01,'FaceAlpha',0.4);
set(gca,'Xscale','log')
title('Control Histogram')
xlabel('Pixel Intensity')
ylabel('Frequency')

%% Control After Induction
% histograms final control distribution

files = dir(fullfile(cpath, tpf_portion));

% compile files into a single table 
for i=a:b
    filepixelmap = readtable(fullfile(cpath, files(i).name));
    [m,n] = size(filepixelmap);
    filepixelmap.Plant = (i+10)*ones(m,1);                                % append plant number as column to table
    if i == a
        totalpixelmap = filepixelmap;
    else    
        totalpixelmap = [totalpixelmap; filepixelmap];
    end
end

% table into arrays
pm = table2array(totalpixelmap);
keys = pm(:,1);         % assigns pixel intensities
values = pm(:,2);       % assigns counts at each intensity

% excludes a plant from histogram data
values(totalpixelmap.Plant == con_exc_plant) = 0;

% compile a list of all unique pixel intensities
ukey = unique(keys);

% sum counts of pixel intensities 
for n = 1:length(ukey)
    totalcounts(n) = sum(values(keys == ukey(n)));
end

% bin pixel intensity counts logarithmically 
for j = 1:(length(logbins)-1)
    loggedcounts(j) = sum(totalcounts((ukey>logbins(j)) & (ukey<logbins(j+1))));
end
% finds the peak of the distribution
[~,I] = max(loggedcounts);
cfpeak = logbins(I);
con_foldchange = cfpeak/cipeak

% add to histogram
hold on
hcf = histogram('BinEdges',logbins,'BinCounts', loggedcounts(1:bins-1),'Normalization','Probability','EdgeAlpha',0.01,'FaceAlpha',0.4);

%% Induced Initial
% histograms initial induced distribution

files = dir(fullfile(ipath, tpi_portion));

for i=a:b
    filepixelmap = readtable(fullfile(ipath, files(i).name));
    [m,n] = size(filepixelmap);

    if i == a
        filepixelmap.Plant = 10*ones(m,1); 
        totalpixelmap = filepixelmap;
    else    
        filepixelmap.Plant = (i-1)*ones(m,1);                                % append plant number as column to table
        totalpixelmap = [totalpixelmap; filepixelmap];
    end
end

% table into arrays
pm = table2array(totalpixelmap);
keys = pm(:,1);         % assigns pixel intensities
values = pm(:,2);       % assigns counts at each intensity

% excludes a plant from histogram data
values(totalpixelmap.Plant == ind_exc_plant) = 0;

% compile a list of all unique pixel intensities
ukey = unique(keys);

% sum counts of pixel intensities 
for n = 1:length(ukey)
    totalcounts(n) = sum(values(keys == ukey(n)));
end

% bin pixel intensity counts logarithmically
for j = 1:(length(logbins)-1)
    loggedcounts(j) = sum(totalcounts((ukey>logbins(j)) & (ukey<logbins(j+1))));
end
% finds the peak of the distribution
[~,I] = max(loggedcounts);
iipeak = logbins(I);

% create histogram
figure()
hci = histogram('BinEdges',logbins,'BinCounts', loggedcounts(1:bins-1),'Normalization','Probability','EdgeAlpha',0.01,'FaceAlpha',0.4);
set(gca,'Xscale','log')
title('Induced Histogram')
xlabel('Pixel Intensity')
ylabel('Frequency')



