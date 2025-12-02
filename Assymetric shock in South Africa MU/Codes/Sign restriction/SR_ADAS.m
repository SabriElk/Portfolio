%% This code is written by Romain Houssa in the framework of the Course,
%% Global Macro-financial Linkages, ESL 2024-2025

% The code does the following things:

%1- estimation of a var(p) model 
%2- identify AD and AS shocks and estimates the associated IRF with sign
%restrictions with F. Canova method of pure sign restriction
%3- recover estimated shocks

clear all; % clear all variables from the workspace
close all; % close all figures
clc;  % Clear Command Window
randn('seed',123456789); % Start with the same random number each time


[A0,B]= xlsread('lab2fruk.xls','fr'); % load Fr and UK data. A contains the data where the first column contains real gdp and the second prices series. B contains text: variables names and dates

% Use first diff to make data statio
y=100*diff(log(A0(1:end,:))); % take the growth rates of series
pmax=4; % maximum # of lags used in the VAR. below we select the optimum one
st=(1970+0/4);
en=(2011+2/4);
timeperQ= st:1/4:en;
[T,n]=size(y); % dimension of the dataset


%-----------------------------------------------------------
%construct the var model with iptimal lag obtained above
%----------------------------------------------------------
p=4;
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

st=(1970+0/4);
en=(2013+2/4);
timeperQ= st:1/4:en;

x0=timeperQ(:,p+1):1/4:timeperQ(:,T);
y0=zeros(1,T-p);


figure; 
subplot(1,2,1);
plot(timeperQ(p+1:T),SH1(:,1),'--');hold on;
plot(timeperQ(p+1:T),SH1(:,2),'LineWidth', 2.5,'Color','r');hold on;
plot(timeperQ(p+1:T),SH1(:,3),'--');grid;line(x0,y0);hold off;
 set(gca,'XLim',[timeperQ(:,p+1) timeperQ(:,T)],'FontSize',15);xlabel('\bf quarters','FontSize',20)

title('Supply Shocks')

subplot(1,2,2);
plot(timeperQ(p+1:T),SH2(:,1),'--');hold on;
plot(timeperQ(p+1:T),SH2(:,2),'LineWidth', 2.5,'Color','r');hold on;
plot(timeperQ(p+1:T),SH2(:,3),'--');grid;line(x0,y0);hold off;
 set(gca,'XLim',[timeperQ(:,p+1) timeperQ(:,T)],'FontSize',15);xlabel('\bf quarters','FontSize',20)

title('Demand Shocks')

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Now plot IFRs
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    sir=sort(IRFid,4); % irf in difference
    sirc=sort(IRFidc,4);% level
    % estract 68 intervals for the responses using identification uncertainty
    mec11=squeeze(sirc(1,1,:,fix(acp*0.50)));
    upc11=squeeze(sirc(1,1,:,fix(acp*0.84)));
    loc11=squeeze(sirc(1,1,:,fix(acp*0.16)));
    
    mec21=squeeze(sirc(2,1,:,fix(acp*0.50)));
    upc21=squeeze(sirc(2,1,:,fix(acp*0.84)));
    loc21=squeeze(sirc(2,1,:,fix(acp*0.16)));
    
  
    mec12=squeeze(sirc(1,2,:,fix(acp*0.50)));
    upc12=squeeze(sirc(1,2,:,fix(acp*0.84)));
    loc12=squeeze(sirc(1,2,:,fix(acp*0.16)));
    
    mec22=squeeze(sirc(2,2,:,fix(acp*0.50)));
    upc22=squeeze(sirc(2,2,:,fix(acp*0.84)));
    loc22=squeeze(sirc(2,2,:,fix(acp*0.16)));
    
    mer=[mec11,mec12,mec21,mec22];
        
     xrom=1:nPeriods; %this is used for Xasis
 
 

  yrom=zeros(1,nPeriods);
        figure;

    subplot(2,2,1),plot([mec11 upc11 loc11]),axis tight;
    legend('Output')
    grid on;
    title('Supply shock','FontWeight','bold','FontSize',10)
    
    subplot(2,2,2),plot([mec21 upc21 loc21]),axis tight;
    legend('Price')
    grid on;
    title('Supply shock','FontWeight','bold','FontSize',10)
    
    % subplot(2,2,3),plot([mec12 upc12 loc12]),axis tight; line(xrom,yrom);hold off;

    subplot(2,2,3),plot([mec12 upc12 loc12]),axis tight;hold off;
    legend('output responses')
    grid on;
    title('Demand shock','FontWeight','bold','FontSize',10)
    
    subplot(2,2,4),plot([mec22 upc22 loc22]),axis tight;
    legend('Price')
    grid on;
    title('Demand shock',...
        'FontWeight','bold','FontSize',10)
    
     xrom=1:nPeriods; %this is used for Xasis
  yrom=zeros(1,nPeriods);

figure;
   
  plot(1:nPeriods,loc11,'--b','LineWidth',2),axis tight;hold on;
    plot(1:nPeriods,upc11,'--b','LineWidth',2),axis tight;hold on;
  plot(1:nPeriods,mec11,'--r','LineWidth',2);hold off
    legend('Output'), grid on;line(xrom,yrom);
    title('Supply shock','FontWeight','bold','FontSize',10)
         set(gca,'XLim',[1 nPeriods],'FontSize',15);xlabel('\bf quarters','FontSize',20)

figure;
 
  plot(1:nPeriods,loc12,'--b','LineWidth',2),axis tight;hold on;
    plot(1:nPeriods,upc12,'--b','LineWidth',2),axis tight;hold on;
  plot(1:nPeriods,mec12,'--r','LineWidth',2);hold off
    legend('Output'), grid on;line(xrom,yrom);
    title('Demand shock','FontWeight','bold','FontSize',10)
         set(gca,'XLim',[1 nPeriods],'FontSize',15);xlabel('\bf quarters','FontSize',20)

figure;
   
  plot(1:nPeriods,loc21,'--b','LineWidth',2),hold on;
    plot(1:nPeriods,upc21,'--b','LineWidth',2),hold on;
  plot(1:nPeriods,mec21,'--r','LineWidth',2);hold off
    legend('Prices'), grid on;line(xrom,yrom);
    title('Supply shock','FontWeight','bold','FontSize',10)
         set(gca,'XLim',[1 nPeriods],'FontSize',15);xlabel('\bf quarters','FontSize',20)

figure;
   
  plot(1:nPeriods,loc22,'--b','LineWidth',2),hold on;
    plot(1:nPeriods,upc22,'--b','LineWidth',2),hold on;
  plot(1:nPeriods,mec22,'--r','LineWidth',2);hold off
    legend('Prices'), grid on;line(xrom,yrom);
    title('Demand shock','FontWeight','bold','FontSize',10)
         set(gca,'XLim',[1 nPeriods],'FontSize',15);xlabel('\bf quarters','FontSize',20)

    