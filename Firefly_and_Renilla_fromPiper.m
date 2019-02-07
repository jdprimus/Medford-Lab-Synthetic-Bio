%% Firefly and Renilla from Piper Output 
% code to normalize and average wells from protoplast experiment
% Jeremy Primus
% March 5, 2018

function sampleTable = Firefly_and_Renilla_fromPiper(Firefly, Renilla, key)
% Firefly - firefly luciferase output matrix from piper
% Renilla - corresponding renilla luciferase output matrix from piper
% key - matrix of sample IDs

FireflyROIsums = ExtractROISumsCol(Firefly);
RenillaROIsums = ExtractROISumsCol(Renilla);
%%
firefly_sums = sum(FireflyROIsums, 1);             % sum columns of firefly (numerically stack image)
renilla_sums = sum(RenillaROIsums, 1);             % sum columns of renilla (numerically stack image)
%a = 7.8014e+09;                             % this is a new standard to normalize output to X promoters - taken from max value of JDP5 renilla in 11/27/17 experiment
%renilla_norm = renilla_sums/a;
normalized = firefly_sums./renilla_sums;     % normalize firefly output to renilla ouput
sampleTable = array2table(normalized);
sampleTable.Properties.VariableNames = key;

% select ROI sums columns from Piper output table
function ROI_sums = ExtractROISumsCol(matrix)
    k = 1;
    j=1;
    while k < size(matrix,2)
        ROI_sums(:,j) = matrix(:,k);
        k = k+8;
        j = j+1;
    end
end
end
