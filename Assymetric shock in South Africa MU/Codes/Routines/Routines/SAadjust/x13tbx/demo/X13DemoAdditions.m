%% X13DemoAdditions

% This demo file demonstrates the use of some auxiliary functions of the X-13
% Toolbox that are not strictly related to the X-13 algorithm itself. In
% particular, the camplet algorithm is used, as well as the addacf, addpcf, and
% addspectrum programs. Notice that these three programs require additional
% Matlab toolboxes.

%#ok<*CHARTEN>

% turn warnings temporarily off
orig_warning_state = warning('off','all');

% line width
lwidth = 78;
% single and double line
sline = repmat('-',1,lwidth+1);
dline = repmat('=',1,lwidth+1);

% display with wrapped lines and leading space
report = @(s) disp(WrapLines(s,lwidth,' '));

% write heading
clc; disp(dline);
report(['DEMONSTRATION OF X-13 TOOLBOX FOR MATLAB : ', ...
    'run camplet and compare with fixedseas and X-13ARIMA-SEATS']);
report(['This script was developed with MATLAB Version ', ...
    '8.3.0.532 (R2014a)']);
disp(sline)
disp(' ');

%% Load US unemployment data from the FRB St.Louis FRED database

disp('Download data from FRED');

% not seasonally adjusted data
ticker = 'UNRATENSA';
url = sprintf(['https://research.stlouisfed.org/fred2/series/%s', ...
      '/downloaddata/%s.csv'], ticker, ticker);
str = urlread(url);
pos    = strfind(str,[char(13),char(10)]);  % ignore first line
format = ['%s','%f'];
cells  = textscan(str(pos(1)+2:end), format, ...
    'delimiter',',', 'treatAsEmpty',{'.'}, 'CollectOutput',1);
fred.ticker = ticker;
fred.dates  = datenum(cells{1},'yyyy-mm-dd');
fred.data   = cells{2};

% seasonally adjusted data by Census Bureau
ticker = 'UNRATE';
url = sprintf(['https://research.stlouisfed.org/fred2/series/%s', ...
      '/downloaddata/%s.csv'], ticker, ticker);
str = urlread(url);
pos    = strfind(str,[char(13),char(10)]);  % ignore first line
format = ['%s','%f'];
cells  = textscan(str(pos(1)+2:end), format, ...
    'delimiter',',', 'treatAsEmpty',{'.'}, 'CollectOutput',1);
fredsa.ticker = ticker;
fredsa.dates  = datenum(cells{1},'yyyy-mm-dd');
fredsa.data   = cells{2};

%% Run camplet, fixedseas, and X-13

disp('Run seasonal adjustments (X-13, camplet, and fixedseas)');

spec = makespec('ADD','CONST','arima','model','(2 1 0)(0 1 1)','X11', ...
    'FIXED','CAMPLET','series','title','Civilian Unemployment Rate');
x = x13(fred.dates,fred.data,spec);
% add official seasonal adjustment to the resul as variable (sao)
x.addvariable('sao',fredsa.dates,fredsa.data,{'sao'},1, ...
    'Census Bureau seasonal adjustment');

%% Add ACF/PACF and spectra of seasonally adjusted series
% This part requires the signal processing and the econometrics toolboxes)

disp('Add ACF, PACF, and spectra');

x.addacf('sa' ,1,'faa','ACF of fixed seasonally adjusted');
x.addpcf('sa' ,1,'fpa','PACF of fixed seasonally adjusted');
x.addacf('csa',1,'caa','ACF of camplet seasonally adjusted');
x.addpcf('csa',1,'cpa','PACF of camplet seasonally adjusted');
x.addacf('d11',1,'xaa','ACF of X-13 seasonally adjusted');
x.addpcf('d11',1,'xpa','PACF of X-13 seasonally adjusted');
x.addspectrum('sa' ,1,'sfa','Spectrum of fixed seasonal adjustment');
x.addspectrum('csa',1,'sca','Spectrum of camplet seasonal adjustment');
x.addspectrum('d11',1,'sxa','Spectrum of X-13 seasonal adjustment');
x.addspectrum('sao',1,'soa','Spectrum of Census seasonal adjustment');

%% Make some graphs

disp('Make some graphs');

fh = figure('Position',[78 88 1295 770]);
plot(fh,x,'faa','fpa','caa','cpa','xaa','xpa');
fh = figure('Position',[440 68 560 802]);
plot(fh,x,'soa','sxa','sca','sfa','rowwise');
fh = figure('Position',[440 68 560 802]);
plot(fh,x,'sao','d11','csa','sa','rowwise');
fh = figure('Position',[31 438 1163 420]);
plot(fh,x,'sao','d11','csa','sa','combined');
