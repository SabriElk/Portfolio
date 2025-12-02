%% DEMO for X13 Toolbox: fixedseas

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
    'run fixedseas and compare with output of X-13ARIMA-SEATS']);
report(['This script was developed with MATLAB Version ', ...
    '8.3.0.532 (R2014a)']);
disp(sline)

%% Preliminaries

% Select type of trend

% LineSmoothing for Matlab up to 2014a
if verLessThan('matlab', '8.4')
    defaultOptions = {'LineSmoothing','on'};
else
    defaultOptions = cell(0);
end

% get correct path
p = fileparts(mfilename('fullpath'));   % directory of this m-file
% if the section is run with Shift-Ctrl-Enter ...
if isempty(p); p = [cd,'\']; end
% location of graphics files for use with X-13-Graph program
grloc = fullfile(p,'graphics\');

% size for figures with subplots
scSize = get(groot,'ScreenSize');
scWidth = scSize(3); scHeight = scSize(4);
sizeFig = @(h,v) round([0.04*scWidth, scHeight*(1-0.08-v/100)-50, ...
    h/100*scWidth, v/100*scHeight]);

size1 = sizeFig(40,80);
size2 = sizeFig(80,45);
size3 = sizeFig(95,45);
size4 = sizeFig(70,70);
size6 = sizeFig(75,68);
size8 = sizeFig(95,65);
size9 = sizeFig(95,90);

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
    'demonstrating fixedseas']);
report(['This script was developed with MATLAB Version ', ...
    '8.3.0.532 (R2014a)']);
disp(sline);

% Here we select the type of the trend: 'ma', 'hp', 'detrend', 'spline'.
% 'polynomial' is possible, too, but does not work well with spr.
method = 'cma';
trans  = 'add';

%% Artificial Data

report(['We start our demonstration of fixedseas with some artificial data. ' ...
    'These data consist of a linear trend, two sinus cycles with ', ...
    'periodicity of 14 and 20, and some noise. We will try to identify ', ...
    'these periods (as if we didn''t know them already) and then filter ', ...
    'them out.']);

% make artificial data with trend, cycles, and noise
nobs = 500;

trend  = 0.025*(1:nobs)' + 5 - 0.00004*(1:nobs)'.^2;
cycle1 = 0.9 * sin((1:nobs)'*(2*pi)/20);
cycle2 = 1.0 * sin((1:nobs)'*(2*pi)/14);
cycle = cycle1 + cycle2;
resid  = 0.6 * randn(nobs,1);
data   = trend + cycle + resid;
% So we know that data has two periods, 14 and 20. We will try to find
% these periods.

t = clock;      % now
m = mod((t(2)-1:-1:t(2)-nobs),12)+1;
y = [0,-cumsum(diff(m)>0)] + t(1);
d = ones(size(y)) * 28;
dates = datenum(y,m,d);
dates = fliplr(dates)';

% We take a look at the data first.
s = x13([dates,data], x13spec);     % import into x13series just to
                                    % plot it nicely
plot(s);
title('Artificial Data');

figure; spr(data,method);
report(['The graph reveals clear spikes at 14 and 20. The spike at ', ...
    '28 is merely an echo of the one at 14, the one at 40 is an echo of ', ...
    'the one at 20, etc. We start by filtering out the shortest cycle ', ...
    '(at 14) first.']);

% filter out period 14
s = x13([dates,data],x13spec('fixedseas','period',14, ...
    'fixedseas','method',method,'fixedseas','transform',trans));
figure; spr(s.sa.sa,method,trans);
report(['The seasonally adjusted series shows no ', ...
    'spike at 14 anymore, but a clear spike at 20 (and possibly also ', ...
    'one are 40, but that is again merely an echo of the 20-cycle). ', ...
    'We now take out period 20 as well.']);

% filter out period 20 as well
s = x13([dates,data],x13spec('fixedseas','period',[14,20], ...
    'fixedseas','method',method,'fixedseas','transform',trans));
figure; spr(s.sa.sa,s.ir.ir,method,trans);
report(['The seasonally adjusted series and the residuals show no ', ...
    'spikes anymore. Possibly some cycles at high frequency (period = 2) ', ...
    'or at very low frequency (high values of period) are marked, but ', ...
    'these are not reliable and can be dismissed.']);

% add true values as variables
s.addvariable('ttr',dates,trend,'ttr',1,'true artificial trend');
s.addvariable('tsf',dates,cycle,'tsf',1,'true artificial cycle');
s.addvariable('tir',dates,resid,'tir',1,'true artificial noise');

% plot true vs estimated values
figure('Position',size1);
ax = subplot(3,1,1); plot(ax,s,'dat','sa','tr','comb');
ax = subplot(3,1,2); plot(ax,s,'sf','tsf','comb');
ax = subplot(3,1,3); plot(ax,s,'ir','tir','comb');

figure('Position',size2);
%
ax = subplot(1,2,1);
scatter(ax,cycle,s.sf.sf,'.b')
axis square; grid on;
xl = xlim; yl = ylim;
xl = [min(xl(1),yl(1)),max(xl(2),yl(2))]; 
xlim([xl(1),xl(2)]); ylim([xl(1),xl(2)]);
c = corrcoef(cycle,s.sf.sf);
title(['\bf',sprintf('seasonal factor (true vs estimated)\ncorrelation: %6.4f', ...
    c(1,2))]);
xlabel('true cycle');
ylabel('estimated seasonal factor');
%
ax = subplot(1,2,2);
scatter(ax, resid,s.ir.ir, '.b')
axis square; grid on;
xl = xlim; yl = ylim;
xl = [min(xl(1),yl(1)),max(xl(2),yl(2))]; 
xlim([xl(1),xl(2)]); ylim([xl(1),xl(2)]);
ok = ~isnan(s.ir.ir);
c = corrcoef(resid(ok),s.ir.ir(ok));
title(['\bf',sprintf('noise (true vs estimated)\ncorrelations: %6.4f', ...
    c(1,2))]);
xlabel('true errors');
ylabel('estimated residuals');

% Finally, we compare the diffrent trend extraction methods.
% Here we do not specity typearg, so we let fixedseas set a default.
% You can, of course, experiment with fifferent settings by adding
% 'fixedseas','typearg',parameter
sma = x13([dates,data],x13spec('fixedseas','period',[14,20], ...
    'fixedseas','method','cma', 'fixedseas','transform',trans, ...
    'series','name','Centered Moving Average'),'quiet');
shp = x13([dates,data],x13spec('fixedseas','period',[14,20], ...
    'fixedseas','method','hp', 'fixedseas','transform',trans, ...
    'series','name','Hodrick-Prescott'),'quiet');
spo = x13([dates,data],x13spec('fixedseas','period',[14,20], ...
    'fixedseas','method','poly', 'fixedseas','transform',trans, ...
    'series','name','Polynomial'),'quiet');
ssp = x13([dates,data],x13spec('fixedseas','period',[14,20], ...
    'fixedseas','method','spline', 'fixedseas','transform',trans, ...
    'series','name','Cubic Spline'),'quiet');
sdt = x13([dates,data],x13spec('fixedseas','period',[14,20], ...
    'fixedseas','method','detrend', 'fixedseas','transform',trans, ...
    'series','name','Detrend'),'quiet');

figure('Position',size2);
ax = subplot(1,2,1);
plot(ax,sma,shp,spo,ssp,sdt,'tr','comb');
ax = subplot(1,2,2);
plot(ax,sma,shp,spo,ssp,sdt,'sa','comb');

%% Real-Life Data

% US. Federal Highway Administration, Vehicle Miles Traveled
% [TRFVOLUSM227NFWA], retrieved from FRED, Federal Reserve Bank of St.
% Louis https://research.stlouisfed.org/fred2/series/TRFVOLUSM227NFWA/,
% December 31, 2014.
load(fullfile(p,'travel'));
name = 'Miles Traveled';

disp(sline)
fprintf(' We apply fixedseas to real data now.\n\n');
report(['Source and description of data: ',travel.source]);

% import into x13series just to plot it nicely
x = x13([travel.dates,travel.data], x13spec);
% A slightly more efficient way to get the data into an x13series object is
% this:
% x = x13series;
% x.addvariable('dat',travel.dates,travel.data,'dat',1);
% This method avoids calling the x13as program (since we don't want to
% compute anything right now).
plot(x);
title(name);

figure; spr(travel.data,method);
report(['The graph reveals a clear spike at 12, and echos at ', ...
    'multiples of 12. These are monthly data with monthly ', ...
    'seasonality, which we filter out now.']);

fspec = x13spec('fixedseas','method',method);
% fspec = x13spec('fixedseas','method','detrend', ...
%     'fixedseas','typearg',datenum(2006,10,31));
x = x13([travel.dates,travel.data],fspec);
figure; spr(x.sa.sa,x.ir.ir,method);
report('The graph reveals no more spikes.');

report(['Let''s compare the result of fixedseas with the output of ', ...
    'X-13ARIMA-SEATS (see X13DemoTravel.m, whose result we replicate here).']);

% newspec = makespec(trans,'X11','DIAG', ...
%     'arima','model','(0 1 [1 4])(2 1 2)', ...
%     'regression','variables', ['(td easter[15] labor[8] thank[1] ', ...
%         'LS1979.May AO1993.May AO1995.Jan)'], ...
%     'regression','save','(td hol ao ls)', ...
%     'x11','seasonalma','s3x5', ...
%     'series','name',name);

xspec = x13spec( ...
    'transform'   ,'function'   ,'none', ...
    'arima'       ,'model'      ,'(0 1 [1 4])(2 1 2)', ...
	'regression'  ,'variables'  ,['(AO1993.May AO1995.Jan LS1979.May ', ...
        'td easter[15] labor[8] thank[1])'], ...
    'regression'  ,'save'       ,'(hol td ao ls)', ...
    'estimate'    ,'save'       ,'(mdl ref rsd rts est lks)', ...
    'check'       ,'save'       ,'(acf ac2 pcf)', ...
    'spectrum'    ,'save'       ,'(sp0 sp1 sp2 spr st0 st1 st2 str)', ...
    'x11'         ,'seasonalma' ,'s3x5', ...
    'x11'         ,'save'       ,'(d8 d10 d11 d12 d13 d16 e2 e3)', ...
    'series'      ,'name'       ,name);

% We instruct the program now to perform the X13ARIMA-SEATS computations
% (as defined in xspec) and the fixedseas computations (as specified in
% fspec).
x = x13([travel.dates,travel.data], x13spec(xspec,fspec));

% plot the result of fixedseas
figure('Position',size4);
ax = subplot(2,2,1); plot(ax,x,'dat','sa','tr','comb');
ax = subplot(2,2,3); plot(ax,x,'ir');
try
    ok = ~isnan(x.ir.ir);
    subplot(2,2,2); parcorr(x.ir.ir(ok),25);
    subplot(2,2,4); autocorr(x.ir.ir(ok),25);
catch
    report('(Sorry. Econometrics Toolbox is not installed.)');
    subplot(2,2,2); spr(x.sa.sa,x.ir.ir,method);
    title('\bfSPR of .sa and .ir');
end
drawnow;
report(['If you have the econometrics toolbox, the graph shows the ', ...
    '(partial) autocorrelation function of the residuals. It is ', ...
    'obvious that significant correlation at period 12 remains, so ', ...
    'fixedseas was unable to remove all the seasonality.']);

% plot the result of X13ARIMA-SEATS
figure('Position',size4);
ax = subplot(2,2,1); plot(ax,x,'dat','d11','d12','comb');
ax = subplot(2,2,3); plot(ax,x,'rsd','d13','comb');
ax = subplot(2,2,2); plot(ax,x,'pcf');
ax = subplot(2,2,4); plot(ax,x,'acf');

% compare the two
figure('Position',size2);
ax = subplot(1,2,1);
plot(ax,x,'d12','tr','comb')
title('\bftrends')
ax = subplot(1,2,2);
plot(ax,x,'e2','sa','comb')
title('\bfseasonally adjusted series')
legend('X-13 (mod seas adj, e2)','fixedseas')
legend('Location','SouthEast')

% plot(x,'d12','tr','comb')

report(['Comparing the seasonally adjusted series (chart on the right) ', ...
    'of fixedseas with X-13ARIMA-SEATS shows that X-13ARIMA-SEATS does ', ...
    'a much better job of removing the seasonal pattern. Its trend ', ...
    '(chart on the left), however, is a little volatile.']);

%% finish up

disp(dline);

% turn warnings on again (or to whatever state they were)
warning(orig_warning_state);
