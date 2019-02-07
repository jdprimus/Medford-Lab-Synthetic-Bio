%% Positive Feedback Histogram
% Jeremy Primus
% March 26, 2018

%% Description 
% Script to automate import of quantified luciferase data and compile
% into pixel intensity histogram

%% Initialization
%clear all;
cpath = 'C:\Users\Medford Lab\Desktop\QuantifiedZFData\MSA107-18\MSA107-18_Control\Result Summary\HistogramData';  % path to control dataset
ipath = 'C:\Users\Medford Lab\Desktop\QuantifiedZFData\MSA107-18\MSA107-18_50\Result Summary\HistogramData';  % path to control dataset

tpi_portion = 'T000*_1_Leaves*';  % designates the timepoint and portion of plant to histogram for the initial histogram
tpf_portion = 'T003*_1_Leaves*';  % designates the timepoint and portion of plant to histogram for the final histogram

% con_exc_plant = 19;             % plant to exclude from control group - set NaN to exclude none
ind_exc_plant = NaN;              % plant to exclude from control group - set NaN to exclude none

bins = 100;                       % number of bins for histograms
logbins = logspace(1.5,6.5,bins); % sets histogram range
loggedcounts = zeros(1,bins);     % initializes bin counts to 0


%% Control Initial
% histograms initial control distribution
files = dir(fullfile(cpath, tpi_portion));  % extract relevant files

% compile files into single table
for i=1:length(files)
    filepixelmap = readtable(fullfile(cpath, files(i).name));
    [m,n] = size(filepixelmap);
    filepixelmap.Plant = (i+10)*ones(m,1);                                % append plant number as column to table
    if i == 1
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
%values(totalpixelmap.Plant == con_exc_plant ) = 0;

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

%% create histogram
figure()
hci = histogram('BinEdges',logbins,'BinCounts', loggedcounts(1:bins-1),'Normalization','Probability','EdgeAlpha',0.01,'FaceAlpha',0.4);
set(gca,'Xscale','log')
title('MSA107-18 Control (T0, T5)')
xlabel('Pixel Intensity')
ylabel('Frequency')

%% Control After Induction
% histograms final control distribution

files = dir(fullfile(cpath, tpf_portion));

% compile files into a single table 
for i=1:length(files)
    filepixelmap = readtable(fullfile(cpath, files(i).name));
    [m,n] = size(filepixelmap);
    filepixelmap.Plant = (i+10)*ones(m,1);                                % append plant number as column to table
    if i == 1
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
%values(totalpixelmap.Plant == con_exc_plant) = 0;

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
%con_foldchange = cfpeak/cipeak

% add to histogram
hold on
hcf = histogram('BinEdges',logbins,'BinCounts', loggedcounts(1:bins-1),'Normalization','Probability','EdgeAlpha',0.01,'FaceAlpha',0.4);
 %% Induced Initial
% histograms initial induced distribution
files = dir(fullfile(ipath, tpi_portion));

for i=1:length(files)
    filepixelmap = readtable(fullfile(ipath, files(i).name));
    [m,n] = size(filepixelmap);

    if i == 1
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
    sumn(n) = totalcounts(n)*ukey(n);
end
avgpxlintinit = sum(sumn)/sum(totalcounts)

% bin pixel intensity counts logarithmically
for j = 1:(length(logbins)-1)
    loggedcounts(j) = sum(totalcounts((ukey>logbins(j)) & (ukey<logbins(j+1))));
end
% finds the peak of the distribution
[~,I] = max(loggedcounts);
iipeak = logbins(I)

% create histogram
figure()
hci = histogram('BinEdges',logbins,'BinCounts', loggedcounts(1:bins-1),'Normalization','Probability','EdgeAlpha',0.01,'FaceAlpha',0.4);
set(gca,'Xscale','log')
title('MSA107-18 50uM T0')
xlabel('Pixel Intensity')
ylabel('Frequency')

%% Induced After Induction
% histograms final control distribution

files = dir(fullfile(ipath, tpf_portion));

% compile files into a single table 
for i=1:length(files)
    filepixelmap = readtable(fullfile(ipath, files(i).name));
    [m,n] = size(filepixelmap);
    filepixelmap.Plant = (i+10)*ones(m,1);                                % append plant number as column to table
    if i == 1
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
values(totalpixelmap.Plant == ind_exc_plant) = 0;

% compile a list of all unique pixel intensities
ukey = unique(keys);

% sum counts of pixel intensities 
for n = 1:length(ukey)
    totalcounts(n) = sum(values(keys == ukey(n)));
    sumn(n) = totalcounts(n)*ukey(n);
end
avgpxlintfin = sum(sumn)/sum(totalcounts)

% bin pixel intensity counts logarithmically 
for j = 1:(length(logbins)-1)
    loggedcounts(j) = sum(totalcounts((ukey>logbins(j)) & (ukey<logbins(j+1))));
end
% finds the peak of the distribution
[~,I] = max(loggedcounts);
ifpeak = logbins(I)
% ind_foldchange = avgpxlintfin/avgpxlintinit

% add to histogram
hold on
hcf = histogram('BinEdges',logbins,'BinCounts', loggedcounts(1:bins-1),'Normalization','Probability','EdgeAlpha',0.01,'FaceAlpha',0.4);



