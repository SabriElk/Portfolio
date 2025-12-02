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

cd("C:\Users\Ribak\Documents\Master_1\Second quadrimestre\Global macro linkages\Project\Codes\Similarities")

%% Initialize matrix of all countries
p_lag={2,4,1,4,2}
countries = {'Mozambique','Zambia','Swaziland','SouthAfrica','Madagascar'}
num_countries = length(countries)
demand_shocks = cell(num_countries, 1);  
supply_shocks = cell(num_countries, 1); 
time_series = cell(num_countries, 1);

%% Create the loop to gather all IRFs of each country

for i = 1:length(countries)

    country = countries{i}
    [A,B]=xlsread("Similiraties_data.xlsx", countries{i}) % Extract data

    y= diff(log(A(:,1:2)))*100 % take the diff log
    
    % setup lag
    p = p_lag{i}

    % set up timeperiod
    dates = datetime(B(2:end,1), 'InputFormat', 'dd/MM/yyyy')
    dates = dates(1:end-p-1)
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

for i = 1:num_countries
    plot(time_series{i}, supply_shocks{i}, 'DisplayName',countries{i}) % add each country in the plot
end
xlabel('Time (Quarters)');
ylabel('Supply Shock');
title('Supply Shocks Across Countries');
legend('show', 'Location', 'best') ;
grid on;
hold off;


%% PLot demand shock accross country

figure;
hold on;

for i = 1:num_countries
    plot(time_series{i}, demand_shocks{i}, 'DisplayName',countries{i}) % add each country in the plot
end
xlabel('Time (Quarters)');
ylabel('Demand Shock');
title('Demand Shocks Across Countries');
legend('show', 'Location', 'best');
grid on;
hold off;

%% Create a correlation matrix
% Initialize matrix

corr_matrix_demand = []
corr_matrix_supply = []

% Create a correlation matrix for supply shock

for i = 1:num_countries
    for j = 1:num_countries

        % create supply shock correlation matrix
        supply_shocks_i= supply_shocks{i}(1:57);
        supply_shocks_j= supply_shocks{j}(1:57);
        corr_matrix_supply(i,j)=corr(supply_shocks_i,supply_shocks_j)

       
    end
end

% Create a correlation matrix for demand shock

for i = 1:num_countries
    for j = 1:num_countries
       
        demand_shocks_i= demand_shocks{i}(1:57);
        demand_shocks_j= demand_shocks{j}(1:57);
        corr_matrix_demand(i,j)=corr(demand_shocks_i,demand_shocks_j)
   
    end
end

% Extract in an excel file

writematrix(corr_matrix_demand,"correlation_matrix_demand_ZR.xlsx")
writematrix(corr_matrix_supply,"correlation_matrix_supply_ZR.xlsx")

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
b1 = b';
I = eye(n*p);
BETA = [b1(:,3:end);I(1:end-n,1:end)]; %here we exclude the constant term

% check whether the VAR is stationary
EI = eig(BETA); % compute the eigen values of BETA
ROOTS = abs(EI); %take the modulus of eigenvalues
ROOTS = sortrows(ROOTS,1); % sort the root in increasing order
display('Roots of the VAR');
ROOTS

%--------------------------
% impulse response analysis
%-------------------------
nPeriods = 144; %# of horizons for irf
XXir1 = zeros(nPeriods,n*p,n); %initialize the matrix where to store the irf

invIminusB = inv(I-BETA);  %(I-B)^-1
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
        XXir1(j,:,jj) = (BETA*XXir1(j-1,:,jj)')';
    end
end

IRFs = XXir1(:,1:n,:); % Keep only current period responses
shocks = inv(C)*res';


end

