% Normalize Treatment group
% Function to normalize and plot three plants in a test group
function normPlot(varargin)   % pass the function three vectors of luminesence data
[~, n] = size(varargin);      % vectors must be the same size
mins = zeros(n,1);
maxes = ones(n,1);
controlAvg = varargin{1}; %(varargin{10}+varargin{11}+varargin{12})/3;
f = ones(n,1);
for i = 1:n                     % get the minimum and maximum value from each vector
     mins(i) = min(varargin{i});    %% switch n to number of control plants to normalize
     maxes(i) = max(varargin{i});   %% control from 0 to 1 and scale test groups to this      
end                                 %% scale

% determine the min and max of the entire data set
     absmin = min(mins);
     absmax = max(maxes);
     figure
     
for i = 1:n    
 % normalize data
      normLum = (varargin{i} - absmin)/(absmax - absmin); 
 % test luminesence divided by control luminesence  
     lum = varargin{i};
     lumdivControl = lum./controlAvg;
     f(i) = subplot(7,1,i);
 % plots
     %plot(varargin{i})      % to plot raw data
     %plot(normLum, 'red')          % to plot normalized data
     plot(lumdivControl, 'red')    % to plot test data divided by control data 
     line([4 4], [0 20])
     hold on
     axis([0,12,0,20])
     linkaxes(f,'xy');
end
