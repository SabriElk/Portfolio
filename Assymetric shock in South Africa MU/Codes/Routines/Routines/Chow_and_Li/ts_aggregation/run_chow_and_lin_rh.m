% How to use Chow&Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Globalisation and macroeconomic policies ECONM823                                                                   %
% 
% Prof. : Romain Houssa                                                   %
%                          %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% GLO
clear all
close all
clc 

[A1,B1]= xlsread('lab2fruk.xls','uk'); % Load quarterly GDP in order to compare the results

Qgdp = A1(:,1)/10^6;

[A,B]= xlsread('lab2fruk.xls','uk_yearly'); % Load yearly gdp

Ygdp = A/10^6;

[A,B]= xlsread('lab2fruk.xls','uk_indicators'); % Load quarterly indicators (use real, SA data!)

Qexp = A(:,1)/10^6;
Qimp = A(:,2)/10^6;

% Chow&Lin inputs
Y = Ygdp; % yearly gdp
x = [Qimp Qexp]; % the quarterly indicators with size(x,1) = 4*size(Y,1)
ta = 1; % pick 1 because gdp is a flow variable
s = 4; % pick 4 because we want annual to quarterly
type = 1; % estimation with weighted least squares. You can try 0 for maximum likelihood (should be almost identical
CLgdp=chowlin(Y,x,ta,s,type);
CLgdp.y %show the results 

%In growth rates
gCLgdp = 100*diff(log(CLgdp.y));
gQgdp = 100*diff(log(Qgdp));

% Compare Chow&Lin with actual quarterly GDP? 
figure;
plot([1970.25:1/4:2013.75],[gCLgdp,gQgdp])
hold on;
legend('Chow&Lin','True GDP growth rate')
title('Comparing True GDP growth rate to Chow&Lin')
grid;


