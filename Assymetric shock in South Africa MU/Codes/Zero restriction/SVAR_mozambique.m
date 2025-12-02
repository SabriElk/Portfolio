%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Students : Sabri El Kaddouri, Ugo Pedini,                               %
%               Jawad  Masuk and Matilde Puce                             %
%                  Group 3
%
%% Code is made to compute                                                %
%  Svar model and IRF for the whole sample and 2 subperiod  using Zero restriction        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 


clear all; % clear all variables from the workspace
close all; % close all figures
clc;  % Clear Command Window

cd("C:\Users\Ribak\Documents\Master_1\SEcond quadri\Global macro linkages\Project\Zero restriction")

%% Whole period response function
clear all

[A,B]=xlsread("Exported_data.xlsx",'Mozambique')
pmax = 4
y = diff(log(A(:,1:2)))*100; % GDP et Prix en taux de croissance
dates = datetime(B(2:end,1), 'InputFormat', 'dd/MM/yyyy')
dates = dates(1:end-1)
timeperQ=dates
if isrow(timeperQ)
    timeperQ = timeperQ';
end

[IRFs_pre, shocks_pre, p_pre] = run_var_analysis(y, 4, timeperQ)


%% For the first subperiod of Mozambique
clear all

[A,B]=xlsread("Exported_data.xlsx",'Mozambique','A1:C41')
pmax = 4
y = diff(log(A(:,1:2)))*100; % GDP et Prix en taux de croissance
dates = datetime(B(2:end,1), 'InputFormat', 'dd/MM/yyyy')
dates = dates(1:end-1)
timeperQ=dates
if isrow(timeperQ)
    timeperQ = timeperQ';
end

[IRFs_pre, shocks_pre, p_pre] = run_var_analysis(y, 4, timeperQ)


%% For the first subperiod of Mozambique
clear all

[A,B]=xlsread("Exported_data.xlsx",'Mozambique','A41:C999')
pmax = 4
y = diff(log(A(:,1:2)))*100; % GDP et Prix en taux de croissance
dates = datetime(B(2:end,1), 'InputFormat', 'dd/MM/yyyy')
dates = dates(1:end-1)
timeperQ=dates
if isrow(timeperQ)
    timeperQ = timeperQ';
end

[IRFs_pre, shocks_pre, p_pre] = run_var_analysis(y, 4, timeperQ)
%%
function [IRFs, shocks, p] = run_var_analysis(y, pmax, timeperQ)
%------------------------------------
% define useful quantities
%-----------------------------------
[T,n] = size(y); % dimension of the dataset

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%Choose the optimal lag level
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
AICval = zeros(pmax,1); %
ilag = 1;

for ilag = 1:pmax
    y1 = y(ilag+1:T,:);  % dependent variable adjusted for the # of lags
    [T1,n] = size(y1); 
    
    %------------------------------------
    % construct the matrix of regressors
    %-----------------------------------
    x = []; %initilize a matrix without specifying the dimension
    % refer to slide in class 13
    i = 1;
    for i = 1:ilag
        x = [x y(ilag-i+1:end-i,:)]; % stacked regressors (lagged y)
    end;

    const = ones(T1,1); % the constant term
    trend = (1:T1)';
    x = [const trend x]; % add the constant term to the list of regressors

    %--------------------
    % ESTIMATE VAR 
    %--------------------
    invxxprim = inv(x'*x); % x'x matrix
    b = invxxprim*(x'*y1); % ols estimator

    BETA = b';
    res = y1 - x*b;   % residuals of the var model
    vmat = (1/T1)*(res'*res);  %vcov of residuals

    %--------------------
    %compute the value of AIC info criteria
    %--------------------
    AICval(ilag,:) = log(det(vmat)) + (2*ilag*(n^2))/T1;
end

[OptAIC, p] = min(AICval);
display('optimal lag is');
p

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
nPeriods = 60; %# of horizons for irf
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
 %-------------------
 %ploting irf
 %-------------------
 xrom=1:nPeriods; %this is used for Xasis
 
 

  yrom=zeros(1,nPeriods);
  figure;plot(xrom,cumsum(XXir1(:,1,1)), 'LineWidth', 2.5,'Color','b');grid;hold on;title('\bf Output responses','FontSize',25);
    plot(xrom,cumsum(XXir1(:,2,1)), 'LineWidth', 2.5,'Color','m');line(xrom,yrom);hold off;
%      set(gca,'XLim',[1 nPeriods],'FontSize',15);legend('\bf Aggregate Supply Shocks','\bf Aggregate Demand Shocks',30);xlabel('\bf quarters','FontSize',20)
       set(gca,'XLim',[1 nPeriods],'FontSize',15);legend('\bf Aggregate Supply Shocks','\bf Aggregate Demand Shocks');xlabel('\bf quarters','FontSize',20)
     
         figure;plot(xrom,cumsum(XXir1(:,1,2)), 'LineWidth', 2.5,'Color','b');grid;hold on;title('\bf Price responses','FontSize',25);
    plot(xrom,cumsum(XXir1(:,2,2)), 'LineWidth', 2.5,'Color','m');line(xrom,yrom);hold off;
    legend('\bf Aggregate Supply Shocks','\bf Aggregate Demand Shocks')
%     legend('\bf Aggregate Supply Shocks','\bf Aggregate Demand Shocks',30)
     set(gca,'XLim',[1 nPeriods],'FontSize',15);xlabel('\bf quarters','FontSize',20)

   
       figure;
    subplot(1,2,1);   plot(xrom,cumsum(XXir1(:,1,1)), 'LineWidth', 2.5,'Color','b');grid;hold on;title('\bf Output responses','FontSize',25);
    plot(xrom,cumsum(XXir1(:,1,2)), 'LineWidth', 2.5,'Color','m');line(xrom,yrom);hold off;
     set(gca,'XLim',[1 nPeriods],'FontSize',15);
     legend('\bf Aggregate Supply Shocks','\bf Aggregate Demand Shocks');xlabel('\bf quarters','FontSize',20)
     
     subplot(1,2,2);plot(xrom,cumsum(XXir1(:,2,1)), 'LineWidth', 2.5,'Color','b');grid;hold on;title('\bf Price responses','FontSize',25);
    plot(xrom,cumsum(XXir1(:,2,2)), 'LineWidth', 2.5,'Color','m');line(xrom,yrom);hold off;
    legend('\bf Aggregate Supply Shocks','\bf Aggregate Demand Shocks')
     set(gca,'XLim',[1 nPeriods],'FontSize',15);xlabel('\bf quarters','FontSize',20)

     
IRFs = XXir1(:,1:n,:); % Keep only current period responses
shocks = inv(C)*res';

%---------------------------------
% Plotting shocks (corrected time indexing)
%---------------------------------
names = {'\bf supply shocks','\bf demand shocks'};
colors = {'b','m'};

min_length = min(length(timeperQ)-p, size(shocks,2));  
valid_time = timeperQ(p+1:p+min_length);               
valid_shocks = shocks(:,1:min_length);                 


figure;
for jj = 1:n
    if jj==1
        subplot(1,2,jj);
        plot(valid_time, valid_shocks(jj,:)', 'LineWidth',2,'Color','b');  
        line(valid_time, zeros(size(valid_time)), 'Color','k');          
        title(names(1,jj),'FontSize',20);
    else
        subplot(1,2,jj);
        plot(valid_time, valid_shocks(jj,:)', 'LineWidth',2,'Color','m');  
        line(valid_time, zeros(size(valid_time)), 'Color','k');            
        title(names(1,jj),'FontSize',20);
    end
    set(gca,'XLim',[valid_time(1) valid_time(end)],'FontSize',15);         
    xlabel('\bf quarters','FontSize',20);

end
end

