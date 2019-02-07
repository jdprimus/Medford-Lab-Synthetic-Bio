%% ONR Report Line Plots
% Jeremy Primus
% Script to automatically import quantified luciferase data and compile
% into plots for ONR report
%
% June 20, 2017

%% Create Table of Control Plants
% set the path to the control data text files
% recognizes '.txt' extension as relevant data
% one .txt file per timepoint, contains data for all plants

clear all;
hours = [0 2.000005 3.666673 5.166675 7.383345 9.783349 11.98335 14.80002 23.86671 26.03338 28.90005 37.06673 48.95008 56.48343 72.41678]; % timepoint of each observation in hours
path = 'C:\Users\Medford Lab\Downloads\OneDrive_1_6-20-2017\MLB101-SC_Control\Result Summary';  % path to control dataset
files = dir(fullfile(path,'*.txt'));
for i=1:length(files)
    data = readtable(fullfile(path, files(i).name));

    if i == 1
        [m,n] = size(data);
        data.Hours = hours(i)*ones(m,1);       % add hours column to table 
        controltable = data;
    else    
        data.Hours = hours(i)*ones(m,1);       % add hours column to table 
        controltable = [controltable; data];
    end
end

%% Create Table of Induced Plants
% set the path to the induced data text files
% recognizes '.txt' extension as relevant data
% one .txt file per timepoint, contains data for all plants

path = 'C:\Users\Medford Lab\Downloads\OneDrive_1_6-20-2017\MLB101-SC_Induced\Result Summary';  % set path to induced dataset
files = dir(fullfile(path,'*.txt'));
for i=1:length(files)
    data = readtable(fullfile(path, files(i).name));
    if i == 1
        [m,n] = size(data);
        data.Hours = hours(i)*ones(m,1);       % add hours column to table 
        inducedtable = data;
    else    
        data.Hours = hours(i)*ones(m,1);       % add hours column to table  
        inducedtable = [inducedtable; data];
    end
end

%% Create Line Plots from Tables
figure()
hold on
title('All Control Plants')
for i = 1:m
    plot(controltable.Hours(controltable.PlantNo_ == i+10),controltable.WholeLucDensity(controltable.PlantNo_ == i+10))
    % calculate fold change for each plant
    foldchangecont(:,i) = controltable.WholeLucDensity(controltable.PlantNo_ == i+10)/controltable.WholeLucDensity(i);
end
%%
figure()
hold on
title('MLB101-12')
for i = 1:(m/2)
    plot(inducedtable.Hours(inducedtable.PlantNo_ == i),inducedtable.WholeLucDensity(inducedtable.PlantNo_ == i))  
    % calculate fold change for each plant
    foldchangeind(:,i) = inducedtable.WholeLucDensity(inducedtable.PlantNo_ == i)/inducedtable.WholeLucDensity(i+1);
end

%%
figure()
hold on
title('MLB101-13 Induced')
for i = (m/2):m
    plot(inducedtable.Hours(inducedtable.PlantNo_ == i),inducedtable.WholeLucDensity(inducedtable.PlantNo_ == i))  
    % calculate fold change for each plant
    foldchangeind(:,i) = inducedtable.WholeLucDensity(inducedtable.PlantNo_ == i)/inducedtable.WholeLucDensity(i+1);
end

%% 
tp = 15;
for j = 1:tp
    conavg(j) = sum(controltable.WholeLucDensity(controltable.Hours == hours(j)))/m;
    indavg(j) = sum(inducedtable.WholeLucDensity(inducedtable.Hours == hours(j)))/m;
end

figure()
plot(hours, conavg, hours, indavg)
legend('Control','Induced')
title('Induced vs. Control Plant Average')

%%
% experiment removing the highest expressing plant from each dataset
inducedtable.WholeLucDensity(inducedtable.PlantNo_ == 8) = 0;
controltable.WholeLucDensity(controltable.PlantNo_ == 19) = 0;
for j = 1:tp
    corrconavg(j) = sum(controltable.WholeLucDensity(controltable.Hours == hours(j)))/(m-1);
    corrindavg(j) = sum(inducedtable.WholeLucDensity(inducedtable.Hours == hours(j)))/(m-1);
end

figure()
plot(hours, corrconavg,hours, corrindavg)
legend('Control','Induced')
title('Induced vs. Control Plant Average: removed top expressing plant from each')

%%
% all fold change
figure()
title('Control Plants Fold Change')
plot(hours,foldchangecont)
figure()
title('Induced Plants Fold Change')
plot(hours,foldchangeind)

% average fold change
foldconavg = sum(foldchangecont, 2)/m;
foldindavg = sum(foldchangeind,2)/m;
figure()
plot(hours, foldconavg, hours, foldindavg)

