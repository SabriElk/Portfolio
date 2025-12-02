%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Student : Sabri El Kaddouri                                             %
%                                                                         %
%% Code is made to compute GDP quarterly                              %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

clear all
clc
close all

%% Import data set 
% Data are CPI, GDP and GDP deflator

cd('C:\Users\Ribak\Documents\Master_1\Second quadrimestre\Global macro linkages\Project\Codes\RawData') % Change in the right directory where data are located

[A,B]= xlsread('Data GDP.xlsx') % Change the name to match your data set

%% Define vectors

GDP = A(1:18,2);
GDP = GDP/1000000; % reshape in millions dollars
Current_accountIndex = A(1:72,3);

disp(length(Current_accountIndex))
disp(length(GDP))

%% Define time index

Time_GDP = 2005:1:2022
Time_Current_accountIndex = 2005:1/4:2022+0.75


figure;plotyy(Time_GDP,GDP, Time_Current_accountIndex, Current_accountIndex)

%% Introduce chowline method

addpath('C:\Users\Ribak\Documents\Master_1\SEcond quadri\Global macro linkages\Project\Routines\Routines\Chow_and_Li\ts_aggregation')


% Setup variable

Y = GDP;
x = [Current_accountIndex];
ta = 1;
s = 4 ;% Annualy to Quaterly
type = 1 ;


Qgdp=chowlin(Y,x,ta,s,type);

chowlinGDP = Qgdp.y

%figure;plotyy(Time_GDP,GDP,Time_Current_accountIndex,Current_accountIndex);hold on;
figure;plotyy(Time_Current_accountIndex,chowlinGDP,Time_GDP,GDP);legend("GDP annualy","GDP quarterly",'Location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seasonaly adjust Quarterly GDP just computed


