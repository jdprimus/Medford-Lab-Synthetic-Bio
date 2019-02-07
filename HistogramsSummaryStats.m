%% Compile all pixel intensity data
% Jeremy Primus
% August 17, 2018

%% Description
% code to compile all pixel intensity data and generate summary statistics
% for use in ABC algorithm

%% Initialization
clear all;
path = 'C:\Users\Medford Lab\Desktop\QuantifiedZFData';
lines = dir(fullfile(path));
lines = lines(~ismember({lines.name},{'.','..'}));
noOfLines = length(lines);
timepoint = {'T000' 'T001' 'T002' 'T003' 'T004' 'T005'};
noOfTimepoints = length(timepoint);
global logbins bins
bins = 100;                       % number of bins for histograms
logbins = logspace(1.5,6.5,bins); % sets histogram range
subplot_index = reshape(1:36, 6, 6)';
byPlant = 0; 
noOfPlants = 6;
f=1;
%%
for i = 1:noOfLines                  % loop over lines
    treatments = dir(fullfile(path, lines(i).name));                               % Extracts the data containing folders for each treatment
    treatments = treatments(arrayfun(@(x) ~strcmp(x.name(1),'.'), treatments));    % Remove hidden files
    figure()
    l = 1;
    for j = 1:length(treatments)     % loop over treatments
        files = dir(fullfile(path, lines(i).name, treatments(j).name, 'Result Summary', 'HistogramData','*.txt'));  % pixel intensity files
        files = files(arrayfun(@(x) ~strcmp(x.name(1),'.'), files));   % Remove hidden files
        % write function to compile histogram from individual file
        for k = 1:noOfTimepoints
            % grab data files from timepoint
            timepointData = files(arrayfun(@(x) contains(x.name, timepoint{k}), files));
            current_path = fullfile(path, lines(i).name, treatments(j).name, 'Result Summary', 'HistogramData');
            % leaves
            timepointData_Leaves = timepointData(arrayfun(@(x) contains(x.name, 'Leaves'), timepointData));
            pixIntCounts = pixelIntensityCounts(timepointData_Leaves, current_path);
            subplot(6,6,subplot_index(l))
            if byPlant == 1
                for n=1:noOfPlants
                    [allPlants_loggedCounts_Leaves(:,i,j,k), logavgperplant(i,j,k,n)] = logBinningforHistograms(pixIntCounts(pixIntCounts.Plant == n,:));
                    plotHistograms(allPlants_loggedCounts_Leaves(:,i,j,k), treatments(j), timepoint{k});
                    hold on
                    f=f+1;
                    summarystatsLine(f) = {lines(i).name};
                    summarystatsTreatment(f) = {treatments(j).name};
                    summarystatsTimepoint(f) = k;
                    summarystatsPlantNo(f) = n;
                    summarystatsLogMean(f) = logavgperplant(i,j,k,n);
                end
            else
                [allPlants_loggedCounts_Leaves(:,i,j,k), logavgpertreat(i,j,k)] = logBinningforHistograms(pixIntCounts);
                plotHistograms(allPlants_loggedCounts_Leaves(:,i,j,k), treatments(j), timepoint{k}, logavgpertreat(i,j,k))
            end
            l = l+1;
            % roots            
%             timepointData_Roots = timepointData(arrayfun(@(x) contains(x.name, 'Roots'), timepointData));
%             allPlants_loggedCounts_Roots = logBinningforHistograms(timepointData_Roots, current_path);
%             % whole 
%             timepointData_Whole = timepointData(arrayfun(@(x) contains(x.name, 'Whole'), timepointData));
%             allPlants_loggedCounts_Whole = logBinningforHistograms(timepointData_Whole, current_path);
            % append all plants to the same histogram
        end
    end
end
% summarystats = table(summarystatsLine', summarystatsTreatment', summarystatsTimepoint', summarystatsPlantNo', summarystatsLogMean');
% summarystats.Properties.VariableNames = {'Line', 'Treatment', 'Timepoint', 'PlantNo','LogMean'};

%% write function to compile data from .txt files to array
function [totalpixelmap] = pixelIntensityCounts(files, currpath)
    % compile files into single table
    for i=1:length(files)
        filepixelmap = readtable(fullfile(currpath, files(i).name));
        [m,~] = size(filepixelmap);
        filepixelmap.Plant = (i)*ones(m,1);                                % append plant number as column to table
        if i == 1
            totalpixelmap = filepixelmap;
        else    
            totalpixelmap = [totalpixelmap; filepixelmap];
        end       
    end
end

function [loggedcounts, logavg] = logBinningforHistograms(totalpixelmap)
global logbins bins
    loggedcounts = zeros(1,bins);     % initializes bin counts to 0
    % table into arrays
    pm = table2array(totalpixelmap);
    keys = pm(:,1);         % assigns pixel intensities
    values = pm(:,2);       % assigns counts at each intensity

    logavg = ((1/sum(values))*sum(log(keys(keys~=0).^values(keys~=0))));
    
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
    
%     % finds the peak of the distribution
    [~,I] = max(loggedcounts);
    cipeak = logbins(I);
end

function plotHistograms(loggedcounts, treatment, timepoint, logavg)
global logbins bins
    % create histogram
    hci = histogram('BinEdges',logbins,'BinCounts', loggedcounts(1:bins-1),'Normalization','Probability','EdgeAlpha',0.01,'FaceAlpha',0.4);
    hold on
    ax = gca;
    set(gca,'Xscale','log')
    title(strcat(treatment.name, '   ', timepoint))
    xlabel('Pixel Intensity')
    ylabel('Frequency')
    line([exp(logavg) exp(logavg)], ax.YLim)
    hold off
end
