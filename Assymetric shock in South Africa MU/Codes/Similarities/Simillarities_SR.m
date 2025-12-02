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
randn('seed',123456789); % Start with the same random number each time

cd("C:\Users\Ribak\Documents\Master_1\Second quadrimestre\Global macro linkages\Project\Codes\Similarities")

%% Initialize matrix of all countries
p_lag={2,4,1,4,2} % apply optimal lag for each country
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
    dates = dates(1:end-1)
    
    timeperQ=dates
    if isrow(timeperQ)
        timeperQ = timeperQ';
    end
    
    [T,n]=size(y)
    
    
    % call the function
    [IRFidc,SH1,SH2, p] = run_SR_analysis(y, p, timeperQ, T );,

    % Store the shock and timeperiod
    demand_shocks{i} = SH2(:,2)'; 
    supply_shocks{i} = SH1(:,2); 
    time_series{i} = timeperQ(1+p:end,1);

end

%% Plot the result for the Supply shock
figure;
hold on;

for i = 1:num_countries
    plot(time_series{i}, supply_shocks{i}, 'DisplayName',countries{i}) % add each country in the plot
end
h_zero = yline(0, 'k-', 'LineWidth', 0.5);
h_zero.Annotation.LegendInformation.IconDisplayStyle = 'off'; % remove the name in legend
xlabel('Time (Quarters)');
ylabel('Supply Shock');
title('Supply Shocks Across Countries');
legend('show', 'Location', 'best') ;
grid on;
hold off;


%% PLot demand shock accross country for the Demand shock

figure;
hold on;

for i = 1:num_countries
    plot(time_series{i}, demand_shocks{i}, 'DisplayName',countries{i}) % add each country in the plot
end
h_zero = yline(0, 'k-', 'LineWidth', 0.5);
h_zero.Annotation.LegendInformation.IconDisplayStyle = 'off'; % remove the name in legend
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
        supply_shocks_i= supply_shocks{i}(1:57)
        supply_shocks_j= supply_shocks{j}(1:57)
        corr_matrix_supply(i,j)=corr(supply_shocks_i,supply_shocks_j)

       
    end
end

% Create a correlation matrix for demand shock

for i = 1:num_countries
    for j = 1:num_countries
       
        demand_shocks_i= demand_shocks{i}(1:57)'
        demand_shocks_j= demand_shocks{j}(1:57)'
        corr_matrix_demand(i,j)=corr(demand_shocks_i,demand_shocks_j)
   
    end
end

% Extract in an excel table
writematrix(corr_matrix_demand,"correlation_matrix_demand_SR.xlsx")
writematrix(corr_matrix_supply,"correlation_matrix_supply_SR.xlsx")

%%
function [IRFidc, SH1, SH2, p] = run_SR_analysis(y, p, timeperQ, T ) 

%-----------------------------------------------------------
%construct the var model with iptimal lag obtained above
%----------------------------------------------------------

y1=y(p+1:T,:);  % dependent variable adjusted for the # of lags that we loose
[T1,n]=size(y1); 

%------------------------------------
% construct the matrix of regressors
%-----------------------------------
x= []; %initilize a matrix without specifying the dimension
i=1;
    for i = 1:p
        x = [x y(p-i+1:end-i,:)]; % stacked regressors (laacped y). the order is first lag first variable first lag second variable secong lag first second lag second ...
    end;

 
const = ones(T1,1); % the constant term
trend=[1:1:max(size(const))]';
x = [const trend x]; % add the constant term to the list of regressors

%--------------------
% ESTIMATE VAR 
%--------------------

invxxprim = inv(x'*x); % x'x matrix
b = invxxprim*(x'*y1); % ols estimator

BETA=b';
res = y1 - x*b;   % residuals of the var model
vmat = (1/T1)*(res'*res);  %vcov of residuals


%define the matrix of coefficients in companion form
b1=b'
I = eye(n*p);
BETA = [b1(:,3:end);I(1:end-n,1:end)] %here we exclude the constant term


%------------------------------
% COMPUTE FIRST CholeskY IRF
%-------------------------------

nPeriods=60; %# of horizons for irf
XXirc    = zeros(nPeriods,n*p,n); %initialize the matrix where to store the cholesky irf

S=chol(vmat,'lower');
jj = 1; % jj variable
 
 for jj=1:2     
 SHOCK = zeros(n,1);
 SHOCK(jj) = 1; %element jj of impulse =1. In this case the 1st element is shocked 
 j=1;
 XXirc(j,1:n*p,jj) = [(S*SHOCK)' zeros(1,n*(p-1))];  
 for j=2:nPeriods;     
     XXirc(j,1:n*p,jj) = (BETA*XXirc(j-1,1:n*p,jj)')'; %for s  
     j=j+1;
 end
 jj=jj+1;
 end
 
 XXirc1=XXirc(:,1:2,:); %now forget about irf of lag variables
 % 1:period; 2:series; 3: shocks 
 % rearrange XXirc1 in the form      % 1:series; 2:shocks; 3: period 
XXirc2    = zeros(n,n,nPeriods); %initialize the matrix 

 h=1;
 for h=1:nPeriods
   XXirc2(:,:,h)  = XXirc1(h,:,:);
   j=j+1;
 end
 
%---------------------------------------------
%Define some useful parameters to be used below
%---------------------------------------------

draw=1; % to control for total number of draws
acp=0; % to control for the accepted draws
IRFid=[]; % initialize the identified IRF
ho = 4;            % number of quarters over which sign restrictions are imposed
Adraws=500;      % maximum number of accpeted draws for the restrictions


while acp<=Adraws
    %---------------------------------------------------------------
    %1. draw a random matrix  and take appropriate tranformations
    %-------------------------------------------------------------
    A=normrnd(0,1,n,n);
    
    [Q R]=qr(A); % Compute qr decomposition of random matrix A
    for ii=1:n;
        if R(ii,ii)<0
           Q(:,ii)=-Q(:,ii); % normalize the sign
        end
    end
    
    %-------------------------------------------
    % 2. Compute candidate IRF implied by draw in Q
    %-------------------------------------------
    
    for j=1:n  %shock
        for i = 1:nPeriods
            XXirc3(:,j,i) = XXirc2(:,:,i)*Q(:,j);
        end
    end

    %-----------------------------------
    % 3. compute cumulated candidate IRF 
    %------------------------------------
    
    XXirc3cum=zeros(n,n,nPeriods);
    for k=1:nPeriods % horizon
        if k==1  % impact period
            XXirc3cum(:,:,1)=XXirc3(:,:,1);
        else
            for c=1:n,
                for r=1:n,
                    XXirc3cum(c,r,k)=XXirc3cum(c,r,k-1)+XXirc3(c,r,k);
                end;
            end;
        end;
    end;
    
    %--------------------------------------------------------------------------------------------------------------
    % 4. check if the draw satisfies the sign restrictions implied by ADAS model are satisfied for responses to shocks at horizons 1 to ho
    %---------------------------------------------------------------------------------------------------------------
    
    z1=[squeeze(XXirc3cum(1,1,1:ho))'>=0 squeeze(XXirc3cum(2,1,1:ho))'<=0 squeeze(XXirc3cum(1,2,1:ho))'>=0 squeeze(XXirc3cum(2,2,1:ho))'>=0 ]
    
    zz1=sum(z1)
%    return
    if zz1==size(z1,2) % check if the restrictions are satisfied
        disp('I found one accepeted draw')
        
        acp=acp+1;
        for j=1:n % for the shocks
            for i = 1:nPeriods
                IRFid(:,j,i,acp) = XXirc2(:,:,i)*Q(:,j); % Now compute the identified IRF for difference series
            end
        end
    end
    
    draw=draw+1;
    
end

disp('# of draws----- # of accepted draws')
[draw-1, acp]

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Now compute identidied cumulated resp i.e. IRF for series in levels
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

IRFidc=zeros(n,n,nPeriods,acp);
for jjj=1:acp % do this for each accepeted draw
    for k=1:nPeriods,
        if k==1,
            IRFidc(:,:,1,jjj)=IRFid(:,:,1,jjj);
        else
            for c=1:size(S,2),
                for r=1:size(S,1),
                    IRFidc(c,r,k,jjj)=IRFidc(c,r,k-1,jjj)+IRFid(c,r,k,jjj);
                end;
            end;
        end;
    end;
end;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Now compute identified shocks for each draw and plot the results
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

for kk=1:acp % do this for each accepeted draw
    sh(:,:,kk)=inv(IRFid(:,:,1,kk))*res'; % compute orthogonal shocks = inv(RQ)*res
end;

for jj=1:size(sh,3),
    mat_cov(:,:,jj)=(cov(sh(:,:,jj)'));
end;
   
s_sh=sort(sh,3);
sh1_50=squeeze(s_sh(1,:,fix(acp*0.50)));
sh1_84=squeeze(s_sh(1,:,fix(acp*0.84)));
sh1_16=squeeze(s_sh(1,:,fix(acp*0.16)));

sh2_50=squeeze(s_sh(2,:,fix(acp*0.50)));
sh2_84=squeeze(s_sh(2,:,fix(acp*0.84)));
sh2_16=squeeze(s_sh(2,:,fix(acp*0.16)));

SH1=[sh1_16' sh1_50' sh1_84'];
SH2=[sh2_16' sh2_50' sh2_84'];

x0=timeperQ(p+1:T)
x0=x0(1):1/4:x0(end);
y0=zeros(size(x0));



end 