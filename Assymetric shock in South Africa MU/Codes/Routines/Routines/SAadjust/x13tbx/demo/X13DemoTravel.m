%% DEMO for X13 Toolbox: Single Series Run

%% Preliminaries

% get correct path
p = fileparts(mfilename('fullpath'));   % directory of this m-file
% if the section is run with Shift-Ctrl-Enter ...
if isempty(p); p = [cd,'\']; end
% location of graphics files for use with X-13-Graph program
grloc = fullfile(p,'graphics\');

% size for figures with subplots
scSize = get(groot,'ScreenSize');   % size of physical monitor in pixels
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
    'run on a single timeseries']);
report(['This script was developed with MATLAB Version ', ...
    '8.3.0.532 (R2014a)']);
disp(sline)

%% Loading Data

% US. Federal Highway Administration, Vehicle Miles Traveled
% [TRFVOLUSM227NFWA], retrieved from FRED, Federal Reserve Bank of St.
% Louis https://research.stlouisfed.org/fred2/series/TRFVOLUSM227NFWA/,
% December 31, 2014.

load(fullfile(p,'travel'));
report(['Source and discription of data: ',travel.source]);
name = 'miles traveled';

% travel   = fetchdata('TRFVOLUSM227NFWA', 'source','fred');
% travelSA = fetchdata('TRFVOLUSM227SFWA', 'source','fred');
% travel   = fetchdata('TRFVOLUSM227NFWA', 'source','fred', ...
%     'from',travelSA.dates(1));

%% Step 1: Quick and Dirty

disp(sline)
fprintf(' Step 1: ''Quick-and-Dirty''\n\n');

report(['We run a seasonal adjustment with the default parameters ', ...
    'and see what is coming out of it.']);

travel1 = x13([travel.dates,travel.data],makespec,'quiet');

disp(travel1.table('transform'));
report('CONCLUSION: The filtering will be additive.');

disp(travel1.table('d8a'));
report('CONCLUSION: The data are clearly seasonal.');

%% Step 2: Identifying Model

disp(sline)
fprintf(' Step 2: Identifying model\n\n');

report(['''',travel1.arima,''' : PICKMDL has not found an ', ...
    'appropriate model. We will allow PICKMDL a wider menu and use ', ...
    'the mixed3.pml file containing 288 models. This will take a while...']);

spec1 = x13spec(makespec,'pickmdl','file','mixed3.pml');
travel1 = x13([travel.dates,travel.data],spec1,'quiet');

report(['''',travel1.arima,''' : PICKMDL has still not found an ', ...
    'appropriate model.']);

report('CONCLUSION: We will use TRAMO instead.');

spec2 = makespec('ADD','TRAMO','X11','AO','LS', 'DIAG', ...
    'series','name',name);
travel2 = x13([travel.dates,travel.data],spec2,'quiet');

report(['TRAMO finds a rather complicated mixed model: ', ...
    travel2.arima, '.']);

disp(travel2.table('regression'));

report(['CONCLUSION: The coefficients appear fine and the roots do ', ...
    'not indicate cancellation of terms. We can provisionally work ', ...
    'with this specification.']);

report(['However, MA(2) and MA(3) are not significant. We remove ', ...
    'them from the regression: (0 1 [1 4])(2 1 2).']);

spec2 = makespec(spec2,'automdl',[],[], ...
    'arima','model','(0 1 [1 4])(2 1 2)');
travel2 = x13([travel.dates,travel.data],spec2,'quiet');

disp(travel2.table('regression'));

report(['NOTE: Two outliers (a level shift and an additive outlier) ', ...
    'have been identified.']);

fh = figure('Position',size2);
plot(fh,travel2,'acf','pcf','comb');

disp(travel2.table('tukey'));
figure('Position',size3);
ax = subplot(1,3,1); plot(ax,travel2,'spr','str','comb');
ax = subplot(1,3,2); plot(ax,travel2,'sp1','st1','comb');
ax = subplot(1,3,3); plot(ax,travel2,'sp2','st2','comb');

disp(travel2.err);

report(['CONCLUSION: The ACF/PACF look good. The spectra also look ', ...
    'clean at the main frequencies, but do contain peaks at trading ', ...
    'day frequencies.']);

%% Step 3: Calendar dummies

disp(sline)
fprintf(' Step 3: Calendar dummies\n\n');

spec3basic = x13spec(spec2,'automdl',[],[],'arima','model',travel2.arima);
spec3test = makespec(spec3basic,'EASTER','TD');

travel3test = x13([travel.dates,travel.data],spec3test,'quiet');

disp(travel3test.table('regr'));
report(['CONCLUSION: Trading day and Easter effects are significant. ', ...
    'We will keep them. Also, now three outliers are identified. ', ...
    'We will fix them as well.']);

spec3 = makespec(spec3basic, 'NO OUTLIERS', ...
    'regression','variables','(td easter[15] LS1979.May AO1993.May AO1995.Jan)', ...
    'regression','save','(td hol ao ls)');
travel3 = x13([travel.dates,travel.data],spec3,'quiet');

disp(travel3.table('tukey'));
figure('Position',size3);
ax = subplot(1,3,1); plot(ax,travel3,'spr','str','comb');
ax = subplot(1,3,2); plot(ax,travel3,'sp1','st1','comb');
ax = subplot(1,3,3); plot(ax,travel3,'sp2','st2','comb');

report('CONCLUSION: The spectra are clean now.');

%% More dummies...

disp(' Trying more dummies...'); disp(' ');

s11 = makespec(spec3,'regression','variables','(labor[1] thank[1])');
s18 = makespec(spec3,'regression','variables','(labor[1] thank[8])');
s81 = makespec(spec3,'regression','variables','(labor[8] thank[1])');
s88 = makespec(spec3,'regression','variables','(labor[8] thank[8])');

t11 =  x13([travel.dates,travel.data],s11,'quiet');
t18 =  x13([travel.dates,travel.data],s18,'quiet');
t81 =  x13([travel.dates,travel.data],s81,'quiet');
t88 =  x13([travel.dates,travel.data],s88,'quiet');

report(['labor[8] and thank[1] or thank[8] are significant. The ', ...
    'labor[8] and thank[1] combination has a higher likelihood, so ', ...
    'we will include these dummies.']);

spec3 = s81;
travel3 = t81;

%% Step 4: Do the Seasonal Filtering

disp(sline)
fprintf(' Step 4: Performing the seasonal filtering\n\n');

spec4 = makespec(spec3,'X11');
travel4 = x13([travel.dates,travel.data],spec4,'quiet');

figure('Position',size4)
ax = subplot(2,2,1);
plot(ax,travel4,'dat','e2','d12','comb');
ax = subplot(2,2,2);
plot(ax,travel4,'d10','bymonth');
ax = subplot(2,2,3);
plot(ax,travel4,'spr','sp1','sp2','comb');
ax = subplot(2,2,4);
plot(ax,travel4,'d13','span','boxplot');

report(['The seasonal factors are rather stable (the graph on the ', ...
    'top right shows little variation). The decomposition (top left ', ...
    'graph) looks reasonable.']);

fh = figure('Position',size6);
seasbreaks(fh,travel4);

report(['A closer inspection into possible seasonal breaks reveals no ', ...
    'major problems either. This graph shows the seasonal factors and ', ...
    'the SI ratios separately for each month. We do see some quantitatively ', ...
    'important shifts, but they are all slow enough so that the X-11 ', ...
    'procedure can deal with it. We do not need to specify seasonal ', ...
    'breaks in the estimation.']);

%% Step 5: Check Stability

disp(sline)
fprintf(' Step 5: Checking stability (this takes a while...)\n\n');

spec5 = makespec(spec4,'SLIDING','HISTORY');
travel5 = x13([travel.dates,travel.data],spec5,'quiet');

% --- sliding span analysis

figure('Position',size4,'Name',[name,': sliding span analysis']);

ax = subplot(2,2,1);
[~,ax] = plot(ax,travel5,'sfs','selection',[0 0 0 0 1]);
title(ax,'\bfmaximum change SA series (sfs)');
ax = subplot(2,2,2);
[~,ax] = plot(ax,travel5,'chs','selection',[0 0 0 0 1]);
title(ax,'\bfmax change seasonal factor (chs)');

ax = subplot(2,2,3);
[~,ax] = plot(ax,travel5,'sfs','selection',[0 0 0 0 1],'span','boxplot');
title(ax,'\bfmaximum change SA series (sfs)');
ax = subplot(2,2,4);
[~,ax] = plot(ax,travel5,'chs','selection',[0 0 0 0 1],'span','boxplot');
title(ax,'\bfmax change seasonal factor (chs)');

report(['CONCLUSION: The sliding span analysis reveals small changes ', ...
    'of the seasonally adjusted series or the seasonal factors. ', ...
    'The maximum revisions are about 2''000, and the level of ', ...
    'the data is between 100''000 and 250''000, so the revisions ', ...
    'amount to about 1%.']);

% --- stability analysis

figure('Position',size4,'Name',[name,': stability analysis']);

ax = subplot(2,2,1);
plot(ax,travel5,'sar')
title(ax,{'\bfmax % change of final vs','concurrent SA series (sar)'});
% % Note: sar = (final./concurrent-1)*100, where
% final = travel4.sae.Final_SA;
% concurrent = travel4.sae.Conc_SA;
% d = travel4.sar.SA_revision-(final./concurrent-1)*100;
% d is equal to zero, except for numerical noise.
ax = subplot(2,2,3);
plot(ax,travel5,'sar','span','boxplot')
title(ax,{'\bfmax % change of final vs','concurrent SA series (sar)'});

ax = subplot(2,2,2);
plot(ax,travel5,'sar','from',datenum(1985,1,1))
title(ax,{'\bf... since 1985'});
ax = subplot(2,2,4);
plot(ax,travel5,'sar','span','boxplot','from',datenum(1985,1,1))
title(ax,{'\bf... since 1985'});

report(['CONCLUSION: The historical analysis also reveals small ', ...
    'changes of the seasonally adjusted series, except in the ', ...
    'beginning of the sample in the late 70s, early 80s.']);

%% Step 6: Adjust Length Of Filter

disp(sline)
fprintf(' Step 6: Adjusting the length of the seasonal filter\n\n');

disp(travel5.table('d9a'));
report(['The X11 procedure selects the length of the filter ', ...
    'according to the global moving seasonality ratio, GMSR. ', ...
    'For a GMSR above 3.5, X11 selects a 3x5 filter, for GMSR below ', ...
    '2.5 it selects a 3x3 filter. Values between 2.5 and 3.5 are in ', ...
    'a grey area, and I don''t know how the filter is selected then.']);

report(['The GMSR for February indicates 3x3 filter for that month. ', ...
    'January, July, and August are in the grey area. However, we can ', ...
    'marginally increase the stability of the filtering by enforcing ', ...
    'a 3x5 filter for all months.']);

spec6 = makespec(spec5,'x11','seasonalma','s3x5');
travel6 = x13([travel.dates,travel.data],spec6, ...
    'graphicsloc',grloc,'quiet');

% --- sliding span and stability analysis

figure('Position',size6,'Name',[name,': sliding span analysis']);

ax = subplot(2,3,1);
[~,ax] = plot(ax,travel5,travel6,'sfs','selection',[0 0 0 0 1],'comb');
title(ax,'\bfmaximum change SA series (sfs)');

% On my computer, the series I'm looking for is called 'Max___DIFF', but on
% others, strangely, it is called ''Max_0x25_DIFF''. To make this computer-
% independent, I look up the fieldnames.
fn5 = fieldnames(travel5.sfs);
fn6 = fieldnames(travel6.sfs);
ax = subplot(2,3,4);
scatter(ax,travel5.sfs.(fn5{end}),travel6.sfs.(fn6{end}),'.');
hold on; plot(xlim,xlim,'k'); grid on;
xlabel(ax,'sfs spec #5');
ylabel(ax,'sfs spec #6');

ax = subplot(2,3,2);
[~,ax] = plot(ax,travel5,travel6,'chs','selection',[0 0 0 0 1],'comb');
title(ax,'\bfmax change seasonal factor (chs)');

fn5 = fieldnames(travel5.chs);
fn6 = fieldnames(travel6.chs);
ax = subplot(2,3,5);
plot(ax,travel5.chs.(fn5{end}),travel6.chs.(fn6{end}),'.');
hold on; plot(xlim,xlim,'k'); grid on;
xlabel(ax,'chs spec #5');
ylabel(ax,'chs spec #6');

ax = subplot(2,3,3);
[~,ax] = plot(ax,travel5,travel6,'sar','comb');
title(ax,{'\bfmax % change of final vs','concurrent SA series (sar)'});

ax = subplot(2,3,6);
plot(ax,travel5.sar.SA_revision,travel6.sar.SA_revision,'.');
hold on; plot(xlim,xlim,'k'); grid on;
xlabel(ax,'sar spec #5');
ylabel(ax,'sar spec #6');
drawnow;

report(['The difference is minimal, but the largest deviations for ', ...
    'spec #5 are made a bit smaller with spec #6.']);

%% Final Step: specification for production

% remove history and sliding spans
spec7 = makespec(spec6, 'history',[],[], 'slidingspans',[],[]);
travel7 = x13([travel.dates,travel.data],spec7,'quiet');

travel7html = x13([travel.dates,travel.data],spec7,'html','quiet');
report('You can view results with web(travel7html.out).')

% report

disp(travel7.spec);
disp(travel7);

disp(travel7.x2d);
disp(travel7.table('tukey'));

figure('Position',size1,'name','X11')
ax = subplot(3,2,1);
plot(ax,travel7,'dat','e2','d12','comb');
ax = subplot(3,2,3);
plot(ax,travel7,'acf','pcf','comb');
ax = subplot(3,2,5);
plot(ax,travel7,'d13','boxplot','span');
ax = subplot(3,2,2);
plot(ax,travel7,'spr','str','comb');
ax = subplot(3,2,4);
plot(ax,travel7,'sp1','st1','comb');
ax = subplot(3,2,6);
plot(ax,travel7,'sp2','st2','comb');

fh = figure('Position',size6);
seasbreaks(fh,travel7);

report('CONCLUSION: The decomposition appears acceptable.');

% % SEATS cannot deal with four lags (or with missing lags for that
% % matter), so we need to specify a shorter model.
%
% spec8 = makespec(spec7,'TRAMO','SEATS');
% travel8 = x13([travel.dates,travel.data],spec8,'quiet');
% 
% figure('Position',size1,'name','SEATS')
% ax = subplot(3,2,1);
% plot(ax,travel8,'dat','s11','s12','comb');
% ax = subplot(3,2,3);
% plot(ax,travel8,'acf','pcf','comb');
% ax = subplot(3,2,5);
% plot(ax,travel8,'s13','boxplot');
% ax = subplot(3,2,2);
% plot(ax,travel8,'spr','str','comb');
% ax = subplot(3,2,4);
% plot(ax,travel8,'s1s','t1s','comb');
% ax = subplot(3,2,6);
% plot(ax,travel8,'s2s','t2s','comb');

disp(dline);
