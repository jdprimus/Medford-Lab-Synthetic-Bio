function [FinalTable] = compileExperimentData(path)
% Jeremy Primus
% June 25, 2018

%% Description 
% Experiment data should first be quantified using Shell3.m script
% Data should be compiled using Post_Analysis.m
% This script compiles all data from a given experiment into one table 
% written for ABC_ZFparforALL

%%  Create master table
lines = dir(path);                                            % identify lines
lines = lines(arrayfun(@(x) ~strcmp(x.name(1),'.'),lines));   % Remove hidden files
noOfLines = length(lines);

for i = 1:noOfLines                  % loop over lines
    treatments = dir(fullfile(path, lines(i).name));                               % Extracts the data containing folders for each treatment
    treatments = treatments(arrayfun(@(x) ~strcmp(x.name(1),'.'), treatments));    % Remove hidden files
    for j = 1:length(treatments)     % loop over treatments
        files = dir(fullfile(path, lines(i).name, treatments(j).name, 'Result Summary', '*.txt'));
        files = files(arrayfun(@(x) ~strcmp(x.name(1),'.'), files));   % Remove hidden files
        for k = 1:length(files)      % loop over timepoints
            TempTable = readtable(fullfile(path, lines(i).name, treatments(j).name, 'Result Summary', files(k).name));
            [noOfPlants, ~] = size(TempTable);
            % initialize empty cells
            line = cell(noOfPlants, 1);
            treatment = cell(noOfPlants, 1);
            timepoint = cell(noOfPlants,1);
            % create vectors of sample properties
            line(:) = {lines(i).name};
            treatment(:) = {treatments(j).name(11:end)};
            timepoint(:) = {files(k).name(4)};
            timepoint = str2double(timepoint);
            % append these properties to table
            TempTable.Line = line;
            TempTable.Treatment = treatment;
            TempTable.Timepoint = timepoint;
            % append new timepoint information to master table
            if i == 1 && j == 1 && k == 1
                FinalTable = TempTable;
            else
                FinalTable = vertcat(FinalTable, TempTable);
            end

        end
    end
end