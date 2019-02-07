function Codon_Optimizer
%% Codon Optimization
% Jeremy Primus
% January 9, 2019
%
% Description:
% This program takes in DNA sequence and outputs a
% codon-optimized version

% Application for species specific codon optimization of a DNA  coding sequence.  
% Default species is Phaeodactylum tricornutum, although DNA sequence can be 
% optimized for other species with codon usage data at https://www.kazusa.or.jp/codon/.

% update January 9
% based on the description of IDT's codon optimization algorithm, https://www.idtdna.com/pages/education/decoded/article/benefits-of-codon-optimization
% codons with a usage frequency of less than 10% are eliminated.  Updating
% so my algorithm does the same.

%% Global variable declaration
global CodonTable

%% default species codon table (in DNA)
% Species:
% Phaeodactylum tricornutum
CodonTable = readtable('PhaeoTriCodonTable.txt');

%% User input to get data
% main ui figure
Mainuifig = uifigure('Visible','off','Position',[100 100 800 300]);     % Create and then hide the UI as it is being constructed.
Mainuifig.Name = 'Codon Optimizer';                                     % Assign the a name to appear in the window title.

% Species select drop-down
SpeciesSelect = uidropdown(Mainuifig);
SpeciesSelect.Items = {'Phaeodactylum tricornutum', 'new species'};
SpeciesSelect.Position = [10 200 200 25];
SpeciesSelectLabel = uilabel(Mainuifig, 'Position', SpeciesSelect.Position + [0 30 0 0], 'Text', 'Select Species');
SpeciesSelect.ValueChangedFcn = @(SpeciesSelect, event) dropdownChanged(SpeciesSelect);

% Sequence entry
EntryField = uieditfield(Mainuifig); 
EntryField.Position = [10 100 780 25];
EntryField.ValueChangedFcn = @(EntryField, event) textChanged(EntryField, SpeciesSelect);
EntryFieldLabel = uilabel(Mainuifig, 'Position', EntryField.Position + [0 30 0 0], 'Text', 'Paste DNA sequence below: (Must be in frame!)');

Mainuifig.Visible = 'on';                                                   % display figure

%% callback functions
% species selected
function dropdownChanged(dd)
    if strcmp(dd.Value, 'new species')
        msgbox('Appropriately formatted codon usage tables can be obtained from https://www.kazusa.or.jp/codon/.  Select "A style like CodonFrequency output in GCG Wisconsin PackageTM" format.  Copy and paste the codon usage table into the text file.  Save the text file (do not change the name!)')
        if ispc 
            !notepad "newspecies.txt"
        elseif ismac
            !textedit "newspecies.txt"
        end
        msgbox('Return to the application window and enter your DNA sequence')
    end
end

% DNA sequence entered
function textChanged(txt, dd)
    if strcmp(dd.Value, 'new species')
        CodonTable = readtable('newspecies.txt');
        CodonTable.Properties.VariableNames = {'AmAcid', 'Codon', 'Number', 'x_1000', 'Fraction'}; 
    elseif strcmp(dd.Value, 'Phaeodactylum tricornutum')
        CodonTable = readtable('PhaeoTriCodonTable.txt');
    end   
    input_seq = upper(txt.Value);
    output_seq = strings;
    if mod(length(input_seq), 3) ~= 0
        errordlg('Number of nucleotides must be a multiple of three.')
    end
    for i=1:(length(input_seq)/3)       
        codon = input_seq(3*i-2:3*i);                   % ID codon 
        idx = strcmp(CodonTable.Codon, codon);      
        AA{i} = CodonTable.AmAcid{idx};            % ID corresponding AA
    synonymous = strcmp(CodonTable.AmAcid, AA{i}); % ID synonymous codons
        if sum(synonymous) > 0
            syn_codons = CodonTable.Codon(synonymous);
            freq = CodonTable.Number(synonymous)/sum(CodonTable.Number(synonymous));    % calculate codon usage frequency
            cumProb = cumsum(freq(freq > 0.1))/max(cumsum(freq(freq > 0.1)));           % cumulative probability, igoring codons w/ frequency less than 0.1, and renormalizing
            draw = rand;                                                                % draw a random number 0 to 1
            tf = draw > cumProb;
            codonidx = sum(tf)+1;
            syn_codons = syn_codons(freq > 0.1);                                        % eliminate synonymous codons with frequency under 0.1
            new_codon = syn_codons{codonidx};                                           % select new codon based on the draw
        end
        output_seq = strcat(output_seq, new_codon);
        idx2 = strcmp(CodonTable.Codon, codon);  
        AA2{i} = CodonTable.AmAcid{idx2};
        codonmatch(i) = strcmp(AA{i}, AA2{i});
    end
    clipboard('copy',output_seq)
    if sum(codonmatch) == length(codonmatch)
        msgbox('Codon optimized sequence copied to clipboard.')
    end
end
end
