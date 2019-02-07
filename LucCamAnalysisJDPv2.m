%% Treatment Average Time Series Plots
% code to compile results from whole plant experiments testing the OHT
% inducible system with and without positive feedback
% Jeremy Primus
% January 16, 2018

%% Identify and Obtain the Data
clear all
% Identify OS - (from Wenlong Xu)
    if ispc 
        Dir1 = 'C:\'; 
    elseif isunix 
        Dir1 = '/'; 
    end

% Identify folder containing all treatments for a given line
    % initialize lineFolder
    lineFolder = uigetdir(Dir1, 'Select the folder containing all treatement folders for a given line'); 
    
    % identify treatments
    Treatments = dir(lineFolder);
    Treatments = Treatments(arrayfun(@(x) ~strcmp(x.name(1),'.'),Treatments)); % Remove hidden files
    treatmentsNo = size(Treatments, 1);     % number of treatments (ie number of distinct folders in original lineFolder)
   
% experimental setup
    % maybe request values here that can be checked later to cofirm all data is
    % there and formatted correctly
        plantsPerPlate = 6;
    % timepoint mapping 
        % here initialize a vector which maps T1, T2, etc. onto actual times
        tMap = [1 17.167 24 41.17 48 65];
    
%% Create a data table unique to each treatment (incuder concentration)
% Treatment loop cycles through all the treatments associated with the selected line    
     for i = 1:treatmentsNo
         % identify group (A or B) and timepoint
         TreatmentABTimePoint = dir(fullfile(Treatments(i).folder, Treatments(i).name, '*_BF'));                    % Extracts the data containing folders for each treatment
         TreatmentABTimePoint = TreatmentABTimePoint(arrayfun(@(x) ~strcmp(x.name(1),'.'),TreatmentABTimePoint));   % Remove hidden files
         treatmentABTimePointNo = size(TreatmentABTimePoint,1);                                                     % identifies the number of data containing folders in each treatment
         
         %% data loop - this loop cycles through all data containing folders
         % associated with a given treatment
         for j = 1:treatmentABTimePointNo
             % identify all the 'TotalDataAnalysis.txt' files
             Datapoint = dir(fullfile(TreatmentABTimePoint(j).folder, TreatmentABTimePoint(j).name,'*TotalDataAnalysis.txt'));
             Datapoint = Datapoint(arrayfun(@(x) ~strcmp(x.name(1),'.'),Datapoint));                                % Remove hidden files
             TempTable = readtable(fullfile(Datapoint.folder, Datapoint.name));                                     % create a temporary table for each .txt file
             totalLucDensityPerPlate(i,j) = sum(TempTable.TotalLuc/sum(TempTable.Area));
             
             %% initialize values for the while loop to assign Plant No.s to tissues
             plantIndex = 1;                                                % beginning value for plant index
             noOfLeafsRoots = size(TempTable,1);                            % number of tissues to assign a Plant No.
             plantNo = zeros(noOfLeafsRoots,1);                             % create a matrix for plant indexes
             totalPlantLucDensity = zeros(plantsPerPlate,1);
             tissueIndex = 2;                                               % starting index for while loop, index 1 will be dealt with outside the loop
            
             if strcmp(TempTable.Tissue(1), 'Leaf ')                        % from what I have seen all text files start with 'Leaf ', LabelNum = 1
                 plantNo(1) = plantIndex;                                   % first tissue belongs to first plant
             else 
                 fprintf('Error!! Data is not in the right format! (Probably) \n');
             end
             %%     assign Plant No. to each tissue
             while (tissueIndex <= noOfLeafsRoots && plantIndex <= plantsPerPlate)
                 while (strcmp(TempTable.Tissue(tissueIndex), 'Leaf '))     % assign consecutive 'Leafs ' to the same plant
                     plantNo(tissueIndex) = plantIndex;
                     if tissueIndex == noOfLeafsRoots
                         break                      
                     end   
                     tissueIndex = tissueIndex + 1;                         % move to next tissue index
                 end

                 while (strcmp(TempTable.Tissue(tissueIndex), 'Root '))     % assign consecutive 'Roots ' to the same plant
                     plantNo(tissueIndex) = plantIndex;
                     if tissueIndex == noOfLeafsRoots
                         break 
                     end   
                     tissueIndex = tissueIndex+1;                           % move to next tissue index
                 end
                 totalPlantLucDensity(plantIndex) = sum(TempTable.TotalLuc(plantNo == plantIndex))/sum(TempTable.Area(plantNo == plantIndex));
                
                 plantIndex = plantIndex + 1;                               % iterate Plant No.
             end
                          
%%          quantification of individual plants
             % to be done: quantify every plant individually
             
              %TotalLuc = sum(TempTable.TotalLuc(TempTable.plantNo == 1))
%               leafTotalLuc
%              leafTotalArea = [];
%              leafLucDensity = []; 
%              rootTotalLuc = [];
%              rootTotalArea = [];
%              rootLucDensity = [];
%              Plant = table(leafTotalLuc, leafTotalArea, leafLucDensity, rootTotalLuc, rootTotalArea, rootLucDensity);

                 

             plantNo = table(plantNo);
             individualPlantDensities(j,:) = totalPlantLucDensity;
             TempTable = [TempTable plantNo];         % add plantNo to TempTable

             %% identify replicate plates in a treatment, label group
             % to be added later: functionality for more plate
             % replicates... A,B,C, etc.
%              h1 = strfind(Datapoint.name, '_A');
%              h2 = strfind(Datapoint.name, '_B');
%              if h1 > 0 
%                  group(j) = {'A'};
%              elseif h2 > 0 
%                  group(j) = {'B'};
%              end
%              

             % identify timepoint
             k = strfind(Datapoint.name, '_T');
             if k > 0
                 timepoint(j) = {Datapoint.name(k+1:k+2)};
             else
                 % error output
                 fprintf('Error!! The folder containing the luciferase quantification files must \n have a designated timepoint in the filename! Ex. "_T1", "_T4" \n');
             end
             
             %% Data Analysis - quantify every plate as a whole
             % to be done later: pre allocate these vectors
             totalLeafArea(j) = sum(TempTable.Area(strcmp(TempTable.Tissue, 'Leaf ')));
             totalLeafLuc(j) = sum(TempTable.TotalLuc(strcmp(TempTable.Tissue, 'Leaf ')));
             leafLucDensity(j) = totalLeafLuc/totalLeafArea;
             totalRootArea(j) = sum(TempTable.Area(strcmp(TempTable.Tissue, 'Root ')));
             totalRootLuc(j) = sum(TempTable.TotalLuc(strcmp(TempTable.Tissue, 'Root ')));
             rootLucDensity(j) = totalRootLuc/totalRootArea;

        end
         
         %% compile a table of important results for all datapoints generated from whole plate assessments 
         totalPlantDensity = (totalLeafLuc + totalRootLuc)./(totalLeafArea + totalRootArea);
         TreatmentTable = table(timepoint', totalLeafArea', totalLeafLuc', leafLucDensity', totalRootArea', totalRootLuc', rootLucDensity', totalPlantDensity', individualPlantDensities );
         TreatmentTable.Properties.VariableNames = {'Timepoint' 'TotalLeafArea' 'TotalLeafLuc' 'LeafLucDensity' 'TotalRootArea' 'TotalRootLuc' 'RootLucDensity', 'TotalPlantDensity', 'IndividualPlantDensities'}; 
         
         % gather information about timepoints and replicates 
         timePoints = unique(TreatmentTable.Timepoint);     % list of unique timepoints in folder titles
         noOfTimePoints = size(timePoints, 1);              % number of unique timepoints in folder titles
         noOfReplicates = 1;         % number of replicates 
         
%%

       %% average fold changes per treatment 
       for ii = 1:noOfReplicates*plantsPerPlate
           % columns
           foldChangesTimeSeries(ii,:) = (TreatmentTable.IndividualPlantDensities(1 +(ii-1)*noOfTimePoints:(noOfTimePoints*ii))/TreatmentTable.IndividualPlantDensities(1+((ii-1)*noOfTimePoints)))
       end
       
       averageFoldChange = mean(foldChangesTimeSeries) % average fold change for a given treatment, all replicate plants and plates
       figure()
       error = std(foldChangesTimeSeries)/2
       errorbar(tMap, averageFoldChange, error)
       title(Treatments(i).name)
%%
%          % combine the results of any replicate plates
%          % to be added later: functionality for more plate
%          % replicates... A,B,C, etc.noOfReplicates*plantsPerPlate
% 
%          
%          % to be added later: code to confirm timepoints are the same
%          % before combining groups
%          combinedGroupsTotalLeafArea = TreatmentTable.TotalLeafArea(strcmp(TreatmentTable.Group, 'A')) + TreatmentTable.TotalLeafArea(strcmp(TreatmentTable.Group, 'B')); 
%          combinedGroupsTotalLeafLuc = TreatmentTable.TotalLeafLuc(strcmp(TreatmentTable.Group, 'A')) + TreatmentTable.TotalLeafLuc(strcmp(TreatmentTable.Group, 'B'));
%          combinedGroupsLeafLucDensity = combinedGroupsTotalLeafLuc./combinedGroupsTotalLeafArea;
%          combinedGroupsTotalRootArea = TreatmentTable.TotalRootArea(strcmp(TreatmentTable.Group, 'A')) + TreatmentTable.TotalRootArea(strcmp(TreatmentTable.Group, 'B'));
%          combinedGroupsTotalRootLuc = TreatmentTable.TotalRootLuc(strcmp(TreatmentTable.Group, 'A')) + TreatmentTable.TotalRootLuc(strcmp(TreatmentTable.Group, 'B'));
%          combinedGroupsRootLucDensity = combinedGroupsTotalRootLuc./combinedGroupsTotalRootArea;
%          combinedGroupsTotalLuc = (combinedGroupsTotalLeafLuc + combinedGroupsTotalRootLuc);
%          combinedGroupsTotalPlantDensity = (combinedGroupsTotalLeafLuc + combinedGroupsTotalRootLuc)./(combinedGroupsTotalLeafArea + combinedGroupsTotalRootArea);
%          
%          CombinedGroupsTreatmentTable = table(timePoints, combinedGroupsTotalLeafArea,  combinedGroupsTotalLeafLuc,  combinedGroupsLeafLucDensity, combinedGroupsTotalRootArea, combinedGroupsTotalRootLuc, combinedGroupsRootLucDensity, combinedGroupsTotalLuc, combinedGroupsTotalPlantDensity);
%          CombinedGroupsTreatmentTable.Properties.VariableNames = {'Timepoint' 'TotalLeafArea' 'TotalLeafLuc' 'LeafLucDensity' 'TotalRootArea' 'TotalRootLuc' 'RootLucDensity', 'TotalLuc','TotalPlantDensity'}; 
%      
%          platesPerTreatment = size(unique(TreatmentTable.Group),1);
%          timePoints = size(unique(TreatmentTable.Timepoint),1);
%%
         % Generate time series plots
         %hold on;
         figure()
         plot(tMap, TreatmentTable.TotalPlantDensity)%, 'Color', a*[1 0 0])
         title(Treatments(i).name);
         %a = a + 0.5;
     end
     
     
% Create a matrix (struct?) for each folder (treatment) - this is the most
% important part.  Identify the right data format to work with
    % define 'NumberofTimepoints' value
    % define a 'PlantsPerPlate' value
    % define a 'PlatesPerTreatment' Value
    % Begin with '_T1'
        % start with '_A_'
            % 1. When encounter 'Leaf' for first time, define Plant #
            % 2. after encounter 'Root', next encounter of 'Leaf', iterate
            % Plant #
            % continue until Plant # = PlantPerPlate
        % next, go to same timepoint, but '_B_'
            % repeat steps 1 and 2
            % continue until Plant # = PlatesPerTreatment*PlantPerPlate

        % Average of all leaves density at timepoint
                % if 'Leaves', add to vector
                % average vector
        % Standard deviation of all leaves density at timepoint
                % standard dev of vector
        % Average roots density at timepoint
                % if 'Roots', add to vector
                % average vector
        % Standard deviation of roots density at timepoint
                % standard dev of vector
    % iterate timepoint
% Plot the average for each treatment across timepoints 1-end
