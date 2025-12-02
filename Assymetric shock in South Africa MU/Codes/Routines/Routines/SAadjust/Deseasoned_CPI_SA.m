%this code is written by Romain Houssa
%course FD SMIDE 2020-2021
clear all; close all; clc;

%%  open data
 [A,B]= xlsread('dartaNSA.xlsx','cpisa'); % load data with date in the
 
%% first column of the numerical output is data

 dates_dateformat=A(:,1); % this needs to be in number format otherwise you will have to transform it first

%% now name the raw series
 
PROD_RH1=A(:,2);

%% now show the path where the toolbox is located if not in the same directory

addpath(genpath('C:\Users\rhoussa\Documents UNamur\TEACHING\2018_2019\Stabiliz\lab\SAadjust\x13tbx'));

%% implement the x13 SA function on your raw data
results1=x13([dates_dateformat,PROD_RH1]);

%% extract the output needed
PROD_RH1sa=results1.d11.d11; 

%% plot now the results together with the raw data
st=(1980+0/4);
en=(2019+3/4);
timeperQ= st:1/4:en;
figure; plot(timeperQ,[PROD_RH1 PROD_RH1sa]); legend('Raw data','SA data'); grid
