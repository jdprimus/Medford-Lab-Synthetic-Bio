function Fluorcam_Analysis_User_App_v1_1
%% Fluorcam User App Version 1.1
% User App to run analysis on fluorcam data
% Version 1.1
% Jeremy Primus
% September 12, 2017

%% Detailed description
    % Application runs principal component analysis on flourcam data to reduce 73 dimensions to 2
    % uses quadratic discriminant analysis with a diagonal covariance matrix
    % estimate (naive Bayes classifiers) to classify the resulting data points
    % as belonging to group 1 or group 2
    % Also uses linear discrimination analysis to create a linear projection of
    % the components and classify the plants

%% Declare variables within the scope of the app
 clear all
 global h tline tline2 tline3 grp k1i k2i k3i k4i k5i score scorecopy
 Table = table;
 Plants = {};
 outlierstf = 0;
 QDAtf = 0;
 KMeanstf = 0;
 adtline = 1;


%% Construct the components of ImportDataFigure

%  Create and then hide the UI as it is being constructed.
uifig = uifigure('Visible','off','Position',[100 100 640 480]);
% Assign the a name to appear in the window title.
uifig.Name = 'Fluorcam Analysis User App';

% Create tabs
tabgrp = uitabgroup(uifig);
tabgrp.Position = [13 15 616 433];
ImportDataTab = uitab(tabgrp, 'Title', 'Import Data');
DataTab = uitab(tabgrp, 'Title', 'Review Data');

% ImportDataTab Components
    % Instructions label
        label = uilabel(ImportDataTab);
        label.Text = {'Data must be stored in Excel spreadsheet in the required PCA data format'; ...
                  'Plants must contain the proper identifying tag anywhere in column 1:'; ...
                  '     Wildtype untreated: ''UT'''; '     Wildtype treated: '' T'''};
        label.Position = [57 317 412 56];

    % Import Data Button
        ImportDataButton = uibutton(ImportDataTab,'Text','Import Data');
        ImportDataButton.Position = [468 275 100 22];
        ImportDataButton.ButtonPushedFcn = @ImportDataButtonPushed;

    % See example button
        ExSprdSheetButton = uibutton(ImportDataTab, 'Text','See Example Spreadsheet');
        ExSprdSheetButton.Position = [275 275 158 22];
        ExSprdSheetButton.ButtonPushedFcn = @SeeexamplespreadsheetButtonPushed;
    
    % Remove outliers Checkbox
        RemoveoutliersCheckBox = uicheckbox(ImportDataTab,'Text','Remove Outliers?', 'Position', [14 53 118 15]);
        RemoveoutliersCheckBox.ValueChangedFcn = @RemoveoutliersCheckBoxValueChanged;

    % Create NoofTransgenicLinesSpinnerLabel
        NoofTransgenicLinesSpinnerLabel = uilabel(ImportDataTab);
        NoofTransgenicLinesSpinnerLabel.HorizontalAlignment = 'right';
        NoofTransgenicLinesSpinnerLabel.Position = [14 216 134 15];
        NoofTransgenicLinesSpinnerLabel.Text = 'No. of Transgenic Lines';

    % Create NoofTransgenicLinesSpinner
        NoofTransgenicLinesSpinner = uispinner(ImportDataTab);
        NoofTransgenicLinesSpinner.Limits = [1 3];
        NoofTransgenicLinesSpinner.ValueDisplayFormat = '%11.1g';
        NoofTransgenicLinesSpinner.Position = [163 212 100 22];
        NoofTransgenicLinesSpinner.Value = 1;
        
    % Clustering algorithm selector
       ClusteringAlgorithmLabel = uilabel(ImportDataTab);
       ClusteringAlgorithmLabel.Text = {'Select a clustering algorithm:'};
       ClusteringAlgorithmLabel.Position = [420 120 412 56];
        
        % KMeans Checkbox
            KMeansCheckBox = uicheckbox(ImportDataTab,'Text','K-Means', 'Position', [450 130 118 15]);
            KMeansCheckBox.ValueChangedFcn = @KMeansCheckBoxValueChanged;
        
        % QDA Checkbox
            QDACheckBox = uicheckbox(ImportDataTab,'Text','QDA', 'Position', [530 130 118 15]);
            QDACheckBox.ValueChangedFcn = @QDACheckBoxValueChanged;
    
    % Create RunPCAButton
        RunAnalysisButton = uibutton(ImportDataTab, 'push');
        RunAnalysisButton.ButtonPushedFcn = @RunAnalysisButtonPushed;
        RunAnalysisButton.Position = [476 46 100 22];
        RunAnalysisButton.Text = 'Run Analysis';
        
% Make the window visible.
uifig.Visible = 'on';
        
%% Construct Components of Linear Discriminant Analysis Figure
% create the figure and components
LDAFig = figure('Visible','off','NumberTitle','off','Position',[400 200 800 600]);
% Assign the a name to appear in the window title.
LDAFig.Name = 'Linear Discriminant Analysis Results';

% Create tabs
LDAtabgrp = uitabgroup(LDAFig);
        AllParametersLDA = uitab(LDAtabgrp, 'Title', 'All Parameters');
        MeasuredParametersLDA = uitab(LDAtabgrp, 'Title', 'Measured Parameters');
        CalculatedParametersLDA = uitab(LDAtabgrp, 'Title', 'Calculated Parameters');
            
%% Construct Components of Principal Component Analysis Figure

% Create figure and tabs
PlotsFig = figure('Name','PCA Plots','NumberTitle','off','Visible','off','Position',[100 100 800 600]);
PlotsTabs = uitabgroup(PlotsFig);
    AllParametersPlot = uitab(PlotsTabs, 'Title', 'All Parameters');
    MeasuredParametersPlot = uitab(PlotsTabs, 'Title', 'Measured Parameters');
    CalculatedParametersPlot = uitab(PlotsTabs, 'Title', 'Calculated Parameters');

% AllParametersTab components
    AllParametersPlantIDButton = uicontrol(AllParametersPlot, 'Style','pushbutton','String','Plant ID','Position',[650 90 100 22],...
        'Callback',@PlantIDButtonPushed);
    AllParametersPlantIDButtonLabel = uicontrol(AllParametersPlot,'HorizontalAlignment','left', 'Style','text','String',...
        {'Click the "Plant ID" button,'; 'aim the crosshairs over data point and click to identify plant'},'Position',[628 120 150 50]);
    AllParametersAxes = axes(AllParametersPlot,'Position',[0.068 0.165 0.680 0.758]);

% MeasuredParametersTab Components
    MeasuredParametersPlantIDButton = uicontrol(MeasuredParametersPlot, 'Style','pushbutton','String','Plant ID','Position',[650 90 100 22],...
        'Callback',@PlantIDButtonPushed);
    MeasuredParametersPlantIDButtonLabel = uicontrol(MeasuredParametersPlot,'HorizontalAlignment','left', 'Style','text','String',...
        {'Click the "Plant ID" button,'; 'aim the crosshairs over data point and click to identify plant'},'Position',[628 120 150 50]);
    MeasuredParametersAxes = axes(MeasuredParametersPlot,'Position',[0.068 0.165 0.680 0.758]);

% CalculatedParametersTab Components
    CalculatedParametersPlantIDButton = uicontrol(CalculatedParametersPlot, 'Style','pushbutton','String','Plant ID','Position',[650 90 100 22],...
        'Callback',@PlantIDButtonPushed);
    CalculatedParametersPlantIDButtonLabel = uicontrol(CalculatedParametersPlot,'HorizontalAlignment','left', 'Style','text','String',...
        {'Click the "Plant ID" button,'; 'aim the crosshairs over data point and click to identify plant'},'Position',[628 120 150 50]);
    CalculatedParametersAxes = axes(CalculatedParametersPlot,'Position',[0.068 0.165 0.680 0.758]);

%% Callback functions
% functions that execute due to user interaction

% Button pushed function: SeeexamplespreadsheetButton
function SeeexamplespreadsheetButtonPushed(~, ~)
    winopen('PCA format.xlsx')
end

% Button pushed function: ImportDataButton
function ImportDataButtonPushed(~, ~)
     [filename, pathname] = uigetfile('*.xlsx','Select Spreadsheet');
     Table = readtable(fullfile(pathname, filename));
     celldata = table2cell(Table);
     celldata = [Table.Properties.VariableNames; celldata];
     UITable = uitable(DataTab, 'Data', celldata,'Position', [0 0 620 400]);
end

% Value changed function: RemoveoutliersCheckBox
function RemoveoutliersCheckBoxValueChanged(~, ~)
    outlierstf = RemoveoutliersCheckBox.Value;
end

% Value changed function: KMeansCheckBoxCheckBox
function KMeansCheckBoxValueChanged(~, ~)
    KMeanstf = KMeansCheckBox.Value;
end

% Value changed function: QDACheckBoxCheckBox
function QDACheckBoxValueChanged(~, ~)
    QDAtf = QDACheckBox.Value;
end

% Button pushed function: RunPCAButton
function RunAnalysisButtonPushed(~, ~)
    % obtain spinner value
    adtline = NoofTransgenicLinesSpinner.Value;
    dlg_title = 'Input';
    if adtline == 1
        prompt = {'Enter transgenic line 1 identifier:'};
        answer = inputdlg(prompt,dlg_title,[1 40]);
        tline = answer;
    elseif adtline == 2
        prompt = {'Enter transgenic line 1 identifier:', ...
                  'Enter transgenic line 2 identifier:'};
        answer = inputdlg(prompt,dlg_title,[1 40]);
        tline = answer(1);
        tline2 = answer(2);
    elseif adtline == 3
        prompt = {'Enter transgenic line 1 identifier:', ...
              'Enter transgenic line 2 identifier:', ...
              'Enter transgenic line 3 identifier:'};
        answer = inputdlg(prompt,dlg_title,[1 40]);
        tline = answer(1);
        tline2 = answer(2);
        tline3 = answer(3);

    end
    %% Extract table data

        % Plant names
        Plants = Table(:,1);
        Plants = table2cell(Plants);
        Plants = cellstr(Plants);

        % get data
        Alldata = Table(:,3:73);
        Alldata = table2array(Alldata);

        % Identify plants
        k1 = strfind(Plants, 'UT');                 % identify untreated WT
        k2 = strfind(Plants, ' T');                 % identify treated WT
        k3 = strfind(Plants, tline);                % identify transgenics

        k1i = find(~cellfun('isempty',k1));         % get indexes of untreated WT
        k2i = find(~cellfun('isempty',k2));         % get indexes of treated WT
        k3i = find(~cellfun('isempty',k3));         % get indexes of transgenics

        % identify, index, and group additional tline if present
        if adtline == 2 || adtline == 3    
            k4 = strfind(Plants, tline2);
            k4i = find(~cellfun('isempty',k4));
            grp(k4i) = 4;
        end
        % identify, index, and group third tline if present
        if adtline == 3   
            k5 = strfind(Plants, tline3);
            k5i = find(~cellfun('isempty',k5));
            grp(k5i) = 5;
        end

        % group test plants
        grp(k1i) = 1;
        grp(k2i) = 2;
        grp(k3i) = 3;

    %% Run Analysis

        % All Parameters
            w = 1./var(Alldata);
            [~, ~, ~, ~, expl, ~] = pca(Alldata,'VariableWeights',w);

            [~, score, ~, ~, ~, ~] = pca(Alldata, 'NumComponents', 2,'VariableWeights',w);

            % calculate t-squared in reduced parameter
            redtsq = mahal(score, score);
            scorecopy = score;

            % Remove outliers if box is checked
            if outlierstf == 1
                score(redtsq > 12.35,:) = NaN;
            end

            % Discriminant analysis of PCA data
            if QDAtf == 1
                [err,fun] = discriminantAnalysis();
                str = {"Principal Component Analysis: Measured Parameters"; "Percent variance explained by principal components 1-5"; ...
                expl(1:5); 'Misclassification error estimate:'; err};
            else
                fun = 0;
                str = {"Principal Component Analysis: All Parameters"; "Percent variance explained by principal components 1-5"; ...
                expl(1:5)};
            end

            XGrid = 0;
            idx2Region = 0; 
            C = 0;

            % KMeans Analysis
            if KMeanstf == 1
                [XGrid, idx2Region, C] = kMeansShades(score);
            end

            % plot the data 
            axes(AllParametersAxes)
            GeneratePlots(fun, score, XGrid, idx2Region, C)
            AllParametersDataOutput = uicontrol(AllParametersPlot,'HorizontalAlignment','left', 'Style','text','String',...
                str, 'Position', [620 340 160 200]); 

            % run LDA
            linearClassifier(Alldata, AllParametersLDA)

        % Measured Parameters
            clear expl score
            w = 1./var(Alldata(:,1:32));
            [~, ~, ~, ~, expl, ~] = pca(Alldata(:,1:32),'VariableWeights',w);

            [~, score, ~, ~, ~, ~] = pca(Alldata(:,1:32), 'NumComponents', 2,'VariableWeights',w);

            % calculate t-squared in reduced parameter
            redtsq = mahal(score, score);
            scorecopy = score;

            % Remove outliers if box is checked
            if outlierstf == 1
                score(redtsq > 12.35,:) = NaN;
            end

            % Disciminant analysis of PCA data
            if QDAtf == 1
                [err,fun] = discriminantAnalysis();
                str = {"Principal Component Analysis: Measured Parameters"; "Percent variance explained by principal components 1-5"; ...
                expl(1:5); 'Misclassification error estimate:'; err};
            else
                fun = 0;
                str = {"Principal Component Analysis: Measured Parameters"; "Percent variance explained by principal components 1-5"; ...
                expl(1:5)};
            end


            XGrid = 0;
            idx2Region = 0; 
            C = 0;

            % KMeans Analysis
            if KMeanstf == 1
                [XGrid, idx2Region, C] = kMeansShades(score);
            end

            % plot the data 
            axes(MeasuredParametersAxes)
            GeneratePlots(fun, score, XGrid, idx2Region, C)
            MeasuredParametersDataOutput = uicontrol(MeasuredParametersPlot,'HorizontalAlignment','left', 'Style','text','String',...
                str, 'Position', [620 340 160 200]); 

            % run LDA
            linearClassifier(Alldata(:,1:32), MeasuredParametersLDA)                

        % Calculated Parameters
            clear expl score
            w = 1./var(Alldata(:,33:71));
            [~, ~, ~, ~, expl, ~] = pca(Alldata(:,33:71),'VariableWeights',w);

            [~, score, ~, ~, ~, ~] = pca(Alldata(:,33:71), 'NumComponents', 2,'VariableWeights',w);

            % calculate t-squared in reduced parameter
            redtsq = mahal(score, score);
            scorecopy = score;

            % Remove outliers if box is checked
            if outlierstf == 1
                score(redtsq > 12.35,:) = NaN;
            end

            % Disciminant analysis of PCA data
            if QDAtf == 1
                [err,fun] = discriminantAnalysis();
                str = {"Principal Component Analysis: Calculated Parameters"; "Percent variance explained by principal components 1-5"; ...
                expl(1:5); 'Misclassification error estimate:'; err};
            else
                fun = 0;
                str = {"Principal Component Analysis: Measured Parameters"; "Percent variance explained by principal components 1-5"; ...
                expl(1:5)};
            end

            XGrid = 0;
            idx2Region = 0; 
            C = 0;

            % KMeans Analysis
            if KMeanstf == 1
                [XGrid, idx2Region, C] = kMeansShades(score);
            end

            % plot the data 
            axes(CalculatedParametersAxes)
            GeneratePlots(fun, score, XGrid, idx2Region, C)
            CalculatedParametersDataOutput = uicontrol(CalculatedParametersPlot,'HorizontalAlignment','left', 'Style','text','String',...
                str, 'Position', [620 340 160 200]); 

            % run LDA
            linearClassifier(Alldata(:,33:71), CalculatedParametersLDA)      

    PlotsFig.Visible = 'on';
    LDAFig.Visible = 'on';
        
end

% Button pushed function: PlantIDButton
function PlantIDButtonPushed(~, ~)
    helpbox = msgbox('Press "Esc" key to exit "Plant ID" function', 'Tip','help','modal');
    waitfor(helpbox)
    gname(Plants,h);
end
    
%% Other functions
% functions not directly called by user interaction

function GeneratePlots(fun, score, XGrid, idx2Region, C)
% GeneratePlots Generates PCA, KMeans, and QDA plots
    % Plot KMeans Region
    if KMeanstf == 1
        gscatter(XGrid(:,1),XGrid(:,2),idx2Region, [1,0.75,1; 1,0.9,0.65; 0.9,0.4,0.55; 0.9,1,0.9;],'..',[],'off');
        hold on;
        Cent = plot(C(:,1), C(:,2), ' kx', 'MarkerSize', 15, 'LineWidth', 3); % mark centroid of each region
        title('PCA with K-Means Clustering');
    end
    
    % Plot PCA datapoints and define legend entries
    if adtline == 2                  % if two tlines are present
        h = gscatter(score(:,1), score(:,2), grp,['g','r','b','k'],'.', [],'on');
        entries = {'WT in H20', 'WT in Salt treated', tline{1}, tline2{1}};
    elseif adtline == 3              % if three tlines present
        h = gscatter(score(:,1), score(:,2), grp,['g','r','b','k','c'],'.', [],'on');
        entries = {'WT in H20', 'WT in Salt treated', tline{1}, tline2{1}, tline3{1}};
    else
        h = gscatter(score(:,1), score(:,2), grp,['g','r','b'],'.', [],'on');
        entries = {'WT in H20', 'WT in Salt treated', tline{1}};
    end
    hold on;
    
    % Plot QDA Region
    if QDAtf == 1
        he = ezplot(fun, [min(score(:,1))-1, max(score(:,1))+1, min(score(:,2))-1, max(score(:,2)+1) ]);
    % set title and legend
        title('PCA with Quadratic Discriminant Classifier');
        legend([h' he],[entries, 'QDA Region'], 'Location', 'best');
    elseif KMeanstf == 1 
        legend([Cent h'], ['K-Means Region Centroid', entries], 'Location', 'best');
    end
    if QDAtf == 1 && KMeanstf == 1
        title('PCA with Quadratic Discriminant Classifier and K-Means Clustering');
        legend([Cent h' he],['K-Means Region Centroid', entries, 'QDA Region'], 'Location', 'best');
    end    
    
    % axis labels
    xlabel 'PC1';
    ylabel 'PC2';
    hold off;
end

function  [err, fun] = discriminantAnalysis()
% runs quadratic discriminat analysis and classifies sample
    training = scorecopy(grp <= 2,:);
    sample = scorecopy(grp > 2,:);
    group = grp(1:size(training,1))';
    [~,err, ~,~,coeff] = classify(sample, training, group,'diagquadratic');
    K = coeff(1,2).const;
    L = coeff(1,2).linear;
    Q = coeff(1,2).quadratic;
    % Function to compute K + L*v + v'*Q*v for multiple vectors
    % v=[x;y]. Accepts x and y as scalars or column vectors.
    fun = @(x,y) K + [x y]*L + sum(([x y]*Q) .* [x y], 2);
end

function  [XGrid, idx2Region, C] = kMeansShades(score)
% kMeansShades 
%   function takes the scores from PCA and assigns regions by kmeans
%   analysis.  Primarily created for the Flourcam PCA user app

    rng(1);             % random number generator
    clusterNo = 4;      % number of clusters to identify
    RepNo = 5;          % number of replicates
    
    % opts = statset('Display','final');
    
    [~,C] = kmeans(score,clusterNo,'Distance','Cityblock','Replicates',RepNo);

    x1 = min(score(:,1))-0.1:0.01:max(score(:,1))+0.1;
    x2 = min(score(:,2))-0.1:0.01:max(score(:,2))+0.1;
    [x1G,x2G] = meshgrid(x1,x2);
    XGrid = [x1G(:),x2G(:)]; % Defines a fine grid on the plot

    idx2Region = kmeans(XGrid,4,'MaxIter',1,'Start',C);     % Assigns each node in the grid to the closest centroid
end

function linearClassifier(Data, where)
% Linear Discriminant Analysis
    training = Data(grp <= 2,:);        % training set: treated and untreated wildtype plants
    sample = Data(grp > 2,:);           % sample set: treated transgenics
    group = grp(1:size(training,1))';   
    % linear classifier call
    [class,~, ~,~,coeff] = classify(sample, training, group,'diaglinear');
    % assigns groups for gscatter call
    traingroups(grp == 1) = 1;
    traingroups(grp == 2) = 2;
    transgroups(class == 1) = 3;
    transgroups(class == 2) = 4;
    finalgroups = [traingroups transgroups];
    
    % coefficients of the linear projection
    K = coeff(1,2).const;
    L = coeff(1,2).linear;

    % project the data to the new axis
    transformData = K + Data*L;
    z = zeros(length(transformData),1);

    % sort by the resulting LDA parameter
    [sortedLDA,Ind] = sort(transformData,'descend');

    % find cutoff value for grouping plants
    a = min(transformData(finalgroups == 3));
    b = max(transformData(finalgroups == 4));
    if a > b
        threshold = min(transformData(finalgroups == 3));   % this is only valid if UT > T
    elseif a < b
        threshold = max(transformData(finalgroups == 3));   % this is only valid if T > UT
    else
        threshold = 0;
    end


    % create the figure and components
    hAxes = axes(where,'NextPlot','add',...   
                 'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
                 'XLim', [min(transformData)-10 max(transformData)+10], ...
                 'Color','none');       
    hAxes.YLim = [0 (0.08*(hAxes.XLim(2)- hAxes.XLim(1)))];
    hAxes.YColor = 'none';
    subplot(2,1,1,hAxes);

    % plot the LDA data
    if threshold ~= 0
        axes(hAxes)
        gscatter(transformData, z, finalgroups,[[0 0.8 0];[1 0 0];[0 1 0.749];[1 0.502 0]],'.', [],'off')
        legend({'WT in H20', 'WT in Salt treated', '"Good" Transgenics', '"Bad" Transgenics'}, 'Location','northoutside');
        xlabel(hAxes, '')

        % add threshold arrow
        arrow = annotation('textarrow','String','Threshold');       % store the arrow information
        arrow.Parent = hAxes;                                       % associate the arrow the the current axes
        arrow.X = [threshold threshold];                            % the location in data units
        arrow.Y = fliplr(get(hAxes, 'YLim'));                       % point the arrow down
        arrow.LineWidth  = 2;                                       % make the arrow bolder for the picture
        arrow.HeadWidth  = 10;
        arrow.HeadLength = 1;
    else 
        axes(hAxes)
        gscatter(transformData, z, finalgroups,[[0 0.8 0];[1 0 0];[1 0.502 0]],'.', [],'off')
        legend({'WT in H20', 'WT in Salt treated', '"Bad" Transgenics'}, 'Location','northoutside');
        xlabel(hAxes, '')
    end

    % Create a color mapping
    % color specification in html code
    clrs(finalgroups == 1) = {'#00CC00'};
    clrs(finalgroups == 2) = {'#FF0000'};
    clrs(finalgroups == 3) = {'#00FFBF'};
    clrs(finalgroups == 4) = {'#FF8000'};
    clrs = clrs';

    % Concatenate html strings:
     outHtml = strcat('<html><table border=0 width=100 bgcolor=', ...
        clrs(Ind), ...                                                  % Choose the appropriate color for each cell
       '"><TR><TD>', ...
       cellfun(@num2str,num2cell(sortedLDA),'UniformOutput',false), ... % Convert num data to cell of chars
       '</TD></TR></body></html>');

    % create and display the UItable
    LDATable = uitable(where,'Data',outHtml,'RowName',Plants(Ind),'ColumnName', 'LDA score', 'Position', [300 5 228 330]);
    LDATable.Position(3) = LDATable.Extent(3);

end
end

