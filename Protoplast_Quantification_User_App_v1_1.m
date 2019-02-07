function Protoplast_Quantification_User_App_v1_1
%% Semi-Automated Protoplast Quantification
% Jeremy Primus
% July 31, 2018

%% Description:
% User application to quantify luminescence data obtained from
% transient/stable protoplast assays in a 96-well plate

%% Global variable declaration
global cal_roi_A1 cal_roi_H1 cal_roi_A12 roi_intensity roi_table current_figure_cont exp_summed_image cont_summed_image bf_image_handle exp_roi_intensity cont_roi_intensity r normalizetf NormalizationCheckBox

%% GUI Interface - Obtain Data
% UI Figure
uifig = uifigure('Visible','off','Position',[100 100 640 480]);     %  Create and then hide the UI as it is being constructed.
uifig.Name = 'Protoplast Quantification User App';                  % Assign the a name to appear in the window title.

% Instructions label
InstructionsLabel = uilabel(uifig);
InstructionsLabel.Text = {'Welcome.  This application allows for semi-automated quantification of luficerase output from';...
          'transient protoplast assays.  In order to use this application, the protoplast data must be collected';...
          'in a 96-well plate.  A brightfield image corresponding to each image stack is necessary for alignment.';...
          'An excel spreadsheet, termed an ROI key, should be provided which identifies the sample in each well.'};
InstructionsLabel.Position = [57 317 600 56];

% Import Experimental Data button
ExperimentalDataButton = uibutton(uifig,'Text','Import Data');
ExperimentalDataButton.Position = [468 275 100 22];
ExperimentalDataButton.ButtonPushedFcn = @importExperimentalDataButtonPushed;

% Normalization Checkbox
NormalizationCheckBox = uicheckbox(uifig,'Text','Yes', 'Position', [14 53 118 15]);
NormalizationCheckBox.ValueChangedFcn = @NormalizationCheckBoxValueChanged;
NormalizationCheckBoxLabel = uilabel(uifig, 'Position', [14 83 400 15]);
NormalizationCheckBoxLabel.Text = 'Do you have a constitutively expressed reporter for normalization?';

% See example button
ExSprdSheetButton = uibutton(uifig, 'Text','See Example Spreadsheet');
ExSprdSheetButton.Position = [275 275 158 22];
ExSprdSheetButton.ButtonPushedFcn = @SeeexamplespreadsheetButtonPushed;

uifig.Visible = 'on';                                                   % display figure
normalizetf = false;

%% Callback functions
% import experimental data
function importExperimentalDataButtonPushed(~,~)
    r = 0;
    exp_datapath = uigetdir('C:\','Select folder containing image stack ');
    [exp_bffile, exp_bfpath] = uigetfile('*.tif','Select corresponding brightfield image');
    [roikey_file, roikey_path] = uigetfile('*.xlsx','Select ROI key corresponding to selected 96-well plate');
    roi_table = readtable(fullfile(roikey_path, roikey_file),'ReadVariableNames',false);
    exp_summed_image = stackImages(exp_datapath);
    alignROI(exp_bffile, exp_bfpath);
    current_figure_cont = false;
end

% Value changed function: NormalizationCheckBox
function NormalizationCheckBoxValueChanged(~, ~)
    normalizetf = NormalizationCheckBox.Value;
end

% open spreadsheet with example ROI layout
function SeeexamplespreadsheetButtonPushed(~, ~)
    winopen('96wellROIformat.xlsx')                      % needs to be in matlab directory
end

% Confirm ROI positions
function ConfirmPositions(~,~)
    r = r+1;
    if current_figure_cont == false
        exp_roi_intensity = getPositionsSumIntensity(exp_summed_image);
        if normalizetf == false
            dataOutput()
        elseif normalizetf == true && r < 2
            pause(2)
            cont_datapath = uigetdir('C:\','Select folder containing the transformation control reporter image stack ');
            [cont_bffile, cont_bfpath] = uigetfile('*.tif','Select corresponding brightfield image');
            h = helpdlg('ROI key should be the same.');
            cont_summed_image = stackImages(cont_datapath);
            alignROI(cont_bffile, cont_bfpath);
        elseif normalizetf == true && r > 1
            current_figure_cont = true;
            cont_roi_intensity = getPositionsSumIntensity(cont_summed_image);
            dataOutput()
        end
    end
end

% function to copy displayed table to clipboard
function copyTable(button_origin,~)
    parentFig = button_origin.Parent;
    copy(parentFig.Children(2).Data)
end

%% Other Functions
% numerically stack image
function summed_image = stackImages(datapath)
    imagelist = dir(fullfile(datapath, '*.tif'));                           % list of images in the stack
    no_of_images = length(imagelist);                                       % number of images in the stack

    % sum the image
        for k = 1:no_of_images
            image_name = imagelist(k).name;                                 % file name of current image                  
            imagepath = fullfile(datapath,image_name);                      % path to current image
            current_image = imread(char(imagepath));                        % selects an image for analysis
            if k == 1
                [m,n] = size(current_image);
                summed_image = uint16(zeros(m,n));
            end
            summed_image = summed_image + current_image;
        end
end

% ROI alignment to brightfield GUI
function alignROI(bffile, bfpath)
    bf_image = imread(fullfile(bfpath, bffile));                                    
    image_window = figure('Position', [100 100 1000 700]);
    % Confirm positions button
    ConfirmPositionsButton = uicontrol(image_window, 'Position', [900 80 100 22],  'String', 'Confirm Postions', 'Callback', @ConfirmPositions);
    imshow(bf_image);
    ConfirmPositionsButton.Visible = 'off';
    bf_image_handle = gca;

    % initial roi positions
    position(:,1,1) = [15.113 119.953 100.422 100.422];
   % position(:,8,1) = [49.8965401785715 869.311941964286 95 95];
   % position(:,1,12) = [1193.65147386278 146.166489550046 95 95];

    % User Calibration
    cal_roi_A1 = imellipse(bf_image_handle, position(:,1,1)');
    setFixedAspectRatioMode(cal_roi_A1,true);
    helpdlg('Move the ROI to well A1 (top left).  Resize if needed.  The size you choose here will be applied to all following ROIs.  Double click on the ROI when complete')
    wait(cal_roi_A1);
    cal_size = getPosition(cal_roi_A1);
    position(:,8,1) = [49.8965401785715 869.311941964286 cal_size(3) cal_size(4)];
    setResizable(cal_roi_A1,false);
    helpdlg('Move the ROI to well H1 (bottom left).  Double click ROI when complete')
    cal_roi_H1 = imellipse(bf_image_handle, position(:,8,1)');
    setResizable(cal_roi_H1,false);
    wait(cal_roi_H1);
    helpdlg('Move the ROI to well A12 (top right).  Double click ROI when complete')
    position(:,1,12) = [1193.65147386278 146.166489550046 cal_size(3) cal_size(4)];
    cal_roi_A12 = imellipse(bf_image_handle, position(:,1,12)');
    setResizable(cal_roi_A12,false);
    wait(cal_roi_A12);
    helpdlg('Make any final adjustments to the three calibration ROIs.  Click "Confirm Positions" when finished.')
    ConfirmPositionsButton.Visible = 'on';
end

% Generate all ROIs and obtain intensities
function roi_intensity =  getPositionsSumIntensity(summed_image)
    % get positions of calibration ROIs
    roi_intensity = zeros(8,12);                                            % matrix to save intensity values
    position(:,1,1) = getPosition(cal_roi_A1);
    position(:,8,1) = getPosition(cal_roi_H1);
    position(:,1,12) = getPosition(cal_roi_A12);
    
    % calculate spacings
    row_slope = (position(:,1,12) - position(:,1,1))/11;
    col_slope = (position(:,8,1) - position(:,1,1))/7;
    
    % create first mask and calculate intensity
    roi_mask(1,1,:,:) = createMask(cal_roi_A1);
    roi_intensity(1,1) = sum(summed_image(roi_mask(1,1,:,:)));
    
    % create ROIs, apply mask and calculate intensities
    for i = 1:8
        for j = 1:11
            position(:,i,j+1) = position(:,i,j) + row_slope;
            roi = imellipse(bf_image_handle, position(:,i,j+1)');
            roi_mask(i,j+1,:,:)= createMask(roi);
            roi_intensity(i,j+1) = sum(summed_image(roi_mask(i,j+1,:,:)));
        end
        if i < 8 
            position(:,i+1,1) = position(:,i,1) + col_slope;
            roi = imellipse(bf_image_handle, position(:,i+1,1)');
            roi_mask(i+1,1,:,:)= createMask(roi);
            roi_intensity(i+1,1) = sum(summed_image(roi_mask(i+1,1,:,:)));
        end
    end
end

% Data output
function dataOutput()
    % identify wells
    idx = ~ismissing(roi_table);
    roikey = table2cell(roi_table);
    % identify samples
    sampleIDs = roikey(idx);
    unique_sampleIDs = unique(sampleIDs,'stable');
    if normalizetf == 1
        roi_intensity = exp_roi_intensity./cont_roi_intensity;
    else
        roi_intensity = exp_roi_intensity;
    end
    % average samples
    for p=1:length(unique_sampleIDs)
        dup_idx = strcmp(unique_sampleIDs(p), roikey);
        avg_intensity(p) = mean(roi_intensity(dup_idx));
        stddev_intensity(p) = std(roi_intensity(dup_idx));
        stderr_intensity(p) = std(roi_intensity(dup_idx))/sqrt(numel(roi_intensity(dup_idx)));
    end
    
    % create a bar graph of the results
     figure
     bar(avg_intensity)
     xticks(1:length(unique_sampleIDs))
     xtickangle(90)
     xticklabels(unique_sampleIDs)
     hold on
     errorbar(avg_intensity, stderr_intensity, '.')
     if normalizetf == true
         title('Normalized Luciferase Output')
     else
         title('Luciferase Output')
     end
     
     % analyzed data summary table
     datasummary_table = table(unique_sampleIDs, avg_intensity', stddev_intensity', stderr_intensity');
     datasummary_table.Properties.VariableNames = {'SampleID', 'AvgIntensity', 'StdDevIntensity', 'StdErrIntensity'};
     datasummary_table_fig = uifigure('Name','Data Summary Table','Position', [100 100 1080 700]);
     UItable_datasummary = uitable(datasummary_table_fig,'Data',datasummary_table, 'Position', [10 10 880 680], 'ColumnEditable',false);
     copy_datasummary_button = uibutton(datasummary_table_fig, 'Text','Copy to Clipboard', 'Position', [900 100 158 22]);
     copy_datasummary_button.ButtonPushedFcn = @copyTable;

     % raw data display
     rawdata_table = table(sampleIDs(:), roi_intensity(idx), 'VariableNames', {'SampleID', 'Intensity'});     
     raw_table_fig = uifigure('Name', 'Raw Data Table','Position', [100 100 1080 700]);
     UItable_rawtable = uitable(raw_table_fig,'Data',rawdata_table, 'Position', [10 10 880 680]);
     copy_rawdata_button = uibutton(raw_table_fig, 'Text','Copy to Clipboard', 'Position', [900 100 158 22]);
     copy_rawdata_button.ButtonPushedFcn = @copyTable;
end
end

