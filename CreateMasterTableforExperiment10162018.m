%% Luciferase Camera Data Analysis
% Jeremy Primus
% October 16, 2018

%% Description 
% Experiment data should first be quantified using Shell3.m script
% Data should be compiled using Post_Analysis.m
% This script is a modified version of 'CreateMasterTableforExperiment.m'
% adjusted for the experimental setup of recent MSA103-13/MSA107-18 liquid
% induction experiment
%% Initialization
clear all;
path = 'C:\Users\Medford Lab\Desktop\ZFLiquidInduction\ImageData'; % The path to experimental data, organized as in this location

%%  Create master table
treatments = dir(path);                                                     % identify treatments
treatments = treatments(arrayfun(@(x) ~strcmp(x.name(1),'.'),treatments));  % Remove hidden files
noOfTreatments = length(treatments);

for i = 1:noOfTreatments                  % loop over treatments 
    files = dir(fullfile(path, treatments(i).name, 'Result Summary', '*.txt'));
    files = files(arrayfun(@(x) ~strcmp(x.name(1),'.'), files));   % Remove hidden files
    for k = 1:length(files)      % loop over timepoints
        TempTable = readtable(fullfile(path, treatments(i).name, 'Result Summary', files(k).name));
        [noOfPlants, ~] = size(TempTable);
        % initialize empty cells
        line = cell(noOfPlants, 1);
        treatment = cell(noOfPlants, 1);
        timepoint = cell(noOfPlants,1);
        % create vectors of sample properties
        treatment(:) = {treatments(i).name};
        line(TempTable.PlantNo_ > 8 & TempTable.PlantNo_ < 17) = {'MSA107-18'};
        line(TempTable.PlantNo_ < 9 & TempTable.PlantNo_ > 0) = {'MSA103-13'};


        timepoint(:) = {files(k).name(4)};
        timepoint = str2double(timepoint);
        % append these properties to table
        TempTable.Line = line;
        TempTable.Treatment = treatment;
        TempTable.Timepoint = timepoint;
        % append new timepoint information to master table
        if i == 1 && k == 1
            FinalTable = TempTable;
        else
            FinalTable = vertcat(FinalTable, TempTable);
        end

    end
end


%% Create Various Plots
tMap = [-1.33 13.75 19.833 44.833 67.833 92];   % mapping timepoints onto hours
uTreatments = unique(FinalTable.Treatment);     % identify treatment groups
uTimepoints = unique(FinalTable.Timepoint);     % identify timepoints
%% Leaf Luc Density: time series plot for each treatment, including all individual plants
% how does ordering work?
for m = 1:length(uTreatments)
    figure()
    for n = 1:noOfPlants/2
        plot(tMap, FinalTable.LeavesLucDensity(FinalTable.PlantNo_== n & strcmp(FinalTable.Treatment, uTreatments(m))))
        hold on
    end
    title(['MSA103-13' uTreatments(m)])
    figure()
    for n = (noOfPlants/2)+1:noOfPlants
        plot(tMap, FinalTable.LeavesLucDensity(FinalTable.PlantNo_== n & strcmp(FinalTable.Treatment, uTreatments(m))))
        hold on
    end
    title(['MSA107-18' uTreatments(m)])
end
%% Whole Plant Luc Density: average of all plants, all treatments
for l = 1:noOfTreatments
    figure()
    for m = 1:length(uTreatments)
        for n = 1:length(uTimepoints)
            tpaverage(m,n) = mean(FinalTable.WholeLucDensity(FinalTable.Timepoint == n-1 & strcmp(FinalTable.Treatment, uTreatments(m)) & strcmp(FinalTable.Line, treatments(l).name)));
            tpstddev(m,n) = std(FinalTable.WholeLucDensity(FinalTable.Timepoint == n-1 & strcmp(FinalTable.Treatment, uTreatments(m)) & strcmp(FinalTable.Line, treatments(l).name)));
        end
    end
     plot(tMap, tpaverage')
     hold on
     title(treatments(l).name)
     legend(uTreatments)
     %axis([0 65 0 7e4])
end

%% Total Whole Luc Density: (total luc of all plants/total area of all plants)



%% Average Luc Leaf Density: average of all plants per treatment
figure()
for m = 1:length(uTreatments)
    for n = 1:length(uTimepoints)
        MSA10313_tpaverage(m,n) = mean(FinalTable.LeavesLucDensity(FinalTable.Timepoint == n-1 & strcmp(FinalTable.Treatment, uTreatments(m)) & strcmp(FinalTable.Line, 'MSA103-13')));
        MSA10313_tpstddev(m,n) = std(FinalTable.LeavesLucDensity(FinalTable.Timepoint == n-1 & strcmp(FinalTable.Treatment, uTreatments(m)) & strcmp(FinalTable.Line, 'MSA103-13')));
        MSA10718_tpaverage(m,n) = mean(FinalTable.LeavesLucDensity(FinalTable.Timepoint == n-1 & strcmp(FinalTable.Treatment, uTreatments(m)) & strcmp(FinalTable.Line, 'MSA107-18')));
        MSA10718_tpstddev(m,n) = std(FinalTable.LeavesLucDensity(FinalTable.Timepoint == n-1 & strcmp(FinalTable.Treatment, uTreatments(m)) & strcmp(FinalTable.Line, 'MSA107-18')));
    end
end
 plot(tMap, MSA10313_tpaverage')
 hold on
 title('MSA103-13')
 legend(uTreatments)
 figure()
 plot(tMap, MSA10718_tpaverage')
 hold on
 title('MSA107-18')
 legend(uTreatments)
%% Average Luc Root Density: average of all plants per treatment
uTimepoints = unique(FinalTable.Timepoint);
for l = 1:noOfTreatments
    figure()
    for m = 1:length(uTreatments)
        for n = 1:length(uTimepoints)
            tpaverage(m,n) = mean(FinalTable.RootsLucDensity(FinalTable.Timepoint == n-1 & strcmp(FinalTable.Treatment, uTreatments(m)) & strcmp(FinalTable.Line, treatments(l).name)));
            tpstddev(m,n) = std(FinalTable.RootsLucDensity(FinalTable.Timepoint == n-1 & strcmp(FinalTable.Treatment, uTreatments(m)) & strcmp(FinalTable.Line, treatments(l).name)));
        end
    end
     plot(tMap, tpaverage')
     hold on
     title(treatments(l).name)
     legend(uTreatments)
end
