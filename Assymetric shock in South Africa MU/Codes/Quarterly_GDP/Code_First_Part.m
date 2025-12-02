%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Students : Sabri El Kaddouri, Ugo Pedini,                               %
%               Jawad  Masuk and Matilde Puce                             %
%                  Group 3
%
%% Code is made to compute GDP quarterly and                              %
%  Compute the seasonnaly adjusted CPI and GDP quaterly computed          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

clear all
clc
close all

%% Import data set 
% Data are CPI, GDP and GDP deflator

cd('C:\Users\Ribak\Documents\Master_1\Second quadrimestre\Global macro linkages\Project\Codes\Quarterly_GDP') % Change in the right directory where data are located

%% Extracting South Africa data, 

outputFilename = 'Exported_data.xlsx' % The final document where data are exported
filename = "Raw_Data.xlsx"
sheetName = matlab.lang.makeValidName('South Africa')

[A,B]=xlsread("Raw_Data.xlsx",'Data_SouthAfrica')

% Exctract data

GDP = A(1:end,2);
CPI = A(1:end,3);
date = A(1:end,1);

%% Seasonnaly adjust GDP and CPI for South Africa
addpath("C:\Users\Ribak\Documents\Master_1\Second quadrimestre\Global macro linkages\Project\Codes\Routines\Routines\SAadjust\x13tbx")
    
results_CPI=x13([date,CPI]) % Using the tool to SA the CPI
results_GDP=x13([date,GDP]) % Using the tool to SA the GDP

% Seasonnally adjust GDP

SA_CPI = results_CPI.d11.d11  % Extract the result
SA_GDP = results_GDP.d11.d11

% Fullfill the excel file with south africa data

date = datetime(date, 'ConvertFrom','excel')
date = datestr(date, 'dd/mm/yyyy'); % Convert number into date

title_date = "date"
writematrix(title_date, outputFilename, 'Sheet',sheetName,'Range','A1')
writematrix(date,outputFilename,'Sheet',sheetName,'Range','A2')

title_CPI = "CPI adjusted for South Africa"
writematrix(title_CPI, outputFilename, 'Sheet',sheetName,'Range','C1')
writematrix(CPI,outputFilename,'Sheet',sheetName,'Range','C2')

title_GDP = 'Quarterly GDP adjusted for South Africa'
writematrix(title_GDP, outputFilename, 'Sheet',sheetName,'Range','B1')
writematrix(GDP, outputFilename, 'Sheet', sheetName,'Range','B2'); 


%% Compute GDP quarterly based on one index for each country that has needed

countries = {'Mozambique','Madagascar','Zambia','Swaziland'} % Apply an index to go into different file of each country

for idx = 1:length(countries)
   
    country = countries{idx} % Choosing the country in the index to work with
    sheetname = sprintf('Data_%s', country) 

    [A,B]=xlsread(filename,sheetname) % extract the data from the excel file

    %%
    % Adjust for each country the size of the date
    
    lastvalidrow = find(~isnan(A(:,2)), 1, 'last') % Apply a filter to not exceed the indent for each country

    GDP = A(1:lastvalidrow,2);
    GDP = GDP/1000000; % reshape in millions dollars
    Quarterly_index = A(:,3);

    disp(length(Quarterly_index));
    disp(length(GDP));

    %% Define time index

    firstdate = datetime(A(1,1), "ConvertFrom","excel")
    yearNumber = year(firstdate)

    Time_GDP = yearNumber:1:2022
    Time_MoneySupply_index = yearNumber:1/4:2022+0.75


    %figure;plotyy(Time_GDP,GDP, Time_MoneySupply_index, moneySupply_index);
    %% 

    addpath('C:\Users\Ribak\Documents\Master_1\Second quadrimestre\Global macro linkages\Project\Codes\Routines\Routines\Chow_and_Li\ts_aggregation')

    % Setup variable

    Y = GDP;
    x = [Quarterly_index];
    ta = 1;
    s = 4 ;% Annualy to Quaterly
    type = 1 ;

    Qgdp=chowlin(Y,x,ta,s,type);
    chowlinGDP = Qgdp.y

    %figure;plotyy(Time_GDP,GDP,Time_Current_accountIndex,Current_accountIndex);hold on;
    %figure;plotyy(Time_MoneySupply_index,chowlinGDP,Time_GDP,GDP);legend("GDP annualy","GDP quarterly",'Location','best');
    
    %% Extract data in

    DateNumber = A(:,1)

    %% Locate where the tool is

    addpath("C:\Users\Ribak\Documents\Master_1\Second quadrimestre\Global macro linkages\Project\Codes\Routines\Routines\SAadjust\x13tbx")

    results_GDP=x13([DateNumber,chowlinGDP]) % Using the tool to Seasonnaly adjust Qgdp
    SA_Quarterly_GDP = results_GDP.d11.d11 % Extract the data

    %%
    % Plot to see the result

    M = 01
    D = 01
    t=datetime(Time_GDP',M,D)

    %% Seasonnally Adjust CPI

    addpath("C:\Users\Ribak\Documents\Master_1\Second quadrimestre\Global macro linkages\Project\Codes\Routines\Routines\SAadjust\x13tbx")
    
    lastDate = find(~isnan(A(:,4)), 1, 'last')
    DateNumber_CPI = A(1:lastDate,1)
    CPI = A(1:lastDate,4) % extract the CPI

    results_CPI=x13([DateNumber_CPI,CPI]) % Using the tool to SA the CPI

    SA_CPI = results_CPI.d11.d11  % Extract the result

    DateNumber_CPI = A(1:end,1)
    DateNumber_CPI = datetime(DateNumber_CPI, 'ConvertFrom','excel')
    Date = datestr(DateNumber_CPI, 'dd/mm/yyyy');
    
    %%
    % Export Excel file for each country and their GDP computed

    sheetName = matlab.lang.makeValidName(country);
    date = B(:,1)
    
    title_date = {sprintf('date')}
    writecell(title_date, outputFilename, 'Sheet',sheetName,'Range','A1')
    writematrix(Date,outputFilename,'Sheet',sheetName,'Range','A2')

    title_CPI = {sprintf('CPI adjusted for %s',country)}
    writecell(title_CPI, outputFilename, 'Sheet',sheetName,'Range','C1')
    writematrix(SA_CPI,outputFilename,'Sheet',sheetName,'Range','C2')

    title_GDP = {sprintf('Quarterly GDP adjusted for %s', country)}
    writecell(title_GDP, outputFilename, 'Sheet',sheetName,'Range','B1')
    writematrix(SA_Quarterly_GDP, outputFilename, 'Sheet', sheetName,'Range','B2');  
end

