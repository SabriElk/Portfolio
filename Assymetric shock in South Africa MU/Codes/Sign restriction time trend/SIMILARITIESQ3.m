%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Students : Sabri El Kaddouri, Ugo Pedini,                               %
%               Jawad  Masuk and Matilde Puce                             %
%                  Group 3
%
%% Code is made to compute                                                %
%  Svar model and shock response for the whole sample and 2 subperiod using sign restriction         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We based our code on code provided by Professor Romain Houssa

clear all; % clear all variables from the workspace
close all; % close all figures
clc;  % Clear Command Window

cd('C:\Users\Ribak\Documents\Master_1\Second quadrimestre\Global macro linkages\Project\Codes\Sign restriction time trend')

%% Initialize matrix of all countries

countries = {'Mozambique','Zambia','Swaziland','SouthAfrica','Madagascar'}
num_countries = length(countries)
demand_shocks = cell(num_countries, 1);  
supply_shocks = cell(num_countries, 1); 
time_series = cell(num_countries, 1);

%% Create the loop to gather all IRFs of each country

for i = 1:length(countries)

    country = countries{i}
    [A,B]=xlsread("SimilaritiesQ3.xlsx", countries{i}) % Extract data

    y= log(A(:,1:2)) % take the log
    
    % setup lag
    pmax = 4 
    p = pmax

    % set up timeperiod
    dates = datetime(B(2:end,1), 'InputFormat', 'dd/MM/yyyy')
    dates = dates(1:end-pmax-1)
    timeperQ=dates
    if isrow(timeperQ)
        timeperQ = timeperQ';
    end
    
    
    % call the function
    [IRFs, shocks, p] = run_var_analysis(y , p,timeperQ );

    % Store the shock and timeperiod
    demand_shocks{i} = shocks(1,:)'; 
    supply_shocks{i} = shocks(2,:)'; 
    time_series{i} = timeperQ(:,1);

end

%% Plot the result
figure;
hold on;
for i = 1:length(countries)
    len = min(length(time_series{i}), length(supply_shocks{i}));
    plot(time_series{i}(1:len), supply_shocks{i}(1:len), 'DisplayName', countries{i});
end
legend show;

xlabel('Time (Quarters)');
ylabel('Supply Shock');
title('Supply Shocks Across Countries');
legend('show', 'Location', 'best');
grid on;
hold off;

%% PLot demand shock accross country

figure;
hold on;
for i = 1:length(countries)
    len = min(length(time_series{i}), length(demand_shocks{i}));
    plot(time_series{i}(1:len), demand_shocks{i}(1:len), 'DisplayName', countries{i});
end
legend show;
xlabel('Time (Quarters)');
ylabel('Demand Shock')
title('Demand Shocks Across Countries');
legend('show', 'Location', 'best') ;
grid on;
hold off;



%% Function
function [IRFs, shocks, p] = run_var_analysis(y, p, timeperQ)
%------------------------------------
% define useful quantities
%-----------------------------------
[T,n] = size(y); % dimension of the dataset

%-----------------------------------------------------------
%% Now construct the var model with optimal lag obtained above
%----------------------------------------------------------

y1 = y(p+1:T,:);  % dependent variable adjusted for the # of lags that we loose
[T1,n] = size(y1); 

%------------------------------------
%% construct the matrix of regressors
%-----------------------------------
x = []; %initilize a matrix without specifying the dimension
for i = 1:p
    x = [x y(p-i+1:end-i,:)]; % stacked regressors
end;

const = ones(T1,1); % the constant term
trend = (1:T1)';
trend2 = trend.^2;
x = [const trend x]; % add the constant term to the list of regressors

%--------------------
%% ESTIMATE VAR 
%--------------------
invxxprim = inv(x'*x); % x'x matrix
b = invxxprim*(x'*y1); % ols estimator

BETA = b';
res = y1 - x*b;   % residuals of the var model
vmat = (1/T1)*(res'*res);  %vcov of residuals

%define the matrix of coefficients in companion form
k = size(y, 2);         % number of variables (e.g. 2)
p = 4;                  % number of lags
Kp = k * p;             % total size of companion matrix

% Extract only the lagged VAR coefficients (b1: [k x kp])
A = BETA(:, end-Kp+1:end);  % ✔ BETA is [k x K], so this gives last kp coefficients
                            % 2×8 matrix if k=2, p=4

% Build the companion matrix B (size: kp x kp)
B = zeros(Kp, Kp);
B(1:k, :) = A;                % first k rows: estimated coefficients
B(k+1:end, 1:end-k) = eye(Kp-k);  % identity block in bottom-left

eigVals = eig(B);             % compute eigenvalues
invIminusB = inv(eye(Kp) - B);  % compute (I - B)^-1 safely

%--------------------------
% impulse response analysis
%-------------------------
nPeriods = 144; %# of horizons for irf
XXir1 = zeros(nPeriods,n*p,n); %initialize the matrix where to store the irf

invIminusB = inv(eye(Kp) - B);  % compute (I - B)^-1 safely
DDprim = invIminusB(1:n,1:n)*vmat*invIminusB(1:n,1:n)'; %DD'
D = chol(DDprim(1:n,1:n),'lower'); % Cholesky decomposition
C = inv(invIminusB(1:n,1:n))*D;

% Initialize and compute IRFs
for jj = 1:n  % For each shock
    SHOCK = zeros(n,1);
    SHOCK(jj) = 1;  % Impulse to variable jj
    
    % response
    impact_response = C*SHOCK;  
    XXir1(1,1:n,jj) = impact_response';  
    XXir1(1,n+1:end,jj) = zeros(1,n*(p-1));  
    
    % Propagate response
    for j = 2:nPeriods
        XXir1(j,:,jj) = (B * XXir1(j-1,:,jj)')';  % use companion matrix
    end
end

IRFs = XXir1(:,1:n,:); % Keep only current period responses
shocks = inv(C)*res';


end

