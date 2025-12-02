%% X13DEMOSEAS --- comparison of X-13, X-12, X-11, and seas

%% load the data and perform various seasonal adjustments

% U.S. Bureau of Labor Statistics, Civilian Unemployment Rate [UNRATENSA]
% retrieved from FRED, Federal Reserve Bank of St. Louis
load unemp

specadd = {'STOCK','ADD','NO OUTLIERS'};     % specs common to all runs

[~,xs] = seas([unemp.dates,unemp.data],12,'add');
x1 = x13([unemp.dates,unemp.data],makespec('X11','FIXED','CAMPLET', ...
    'series','title','simplified X-11',specadd{:}), 'x-11','-n');
x2 = x13([unemp.dates,unemp.data], ...
    makespec('DEFAULT','series','title','X-12',specadd{:}), 'x-12');
x3 = x13([unemp.dates,unemp.data], ...
    makespec('X','pickmdl','file','pure4.pml', ...
    'series','title','X-13ARIMA',specadd{:}));
x3s = x13([unemp.dates,unemp.data], ...
    makespec('S','SEATS','series','title','TRAMO-SEATS',specadd{:}));

%% look at the unfiltered data

plot(x3);

%% results of different algorithms

figure('Position',[315 65 667 800]);
ticks = {'multdateticks',10};

ah = subplot(6,3,1); plot(ah,x1,'tr','sa','combined',ticks{:});
ylim([0,12]); title('\rm tr and sa'); ylabel('fixedseas');
ah = subplot(6,3,2); plot(ah,x1,'sf','ir','combined',ticks{:});
title('\rm sf and ir');
ah = subplot(6,3,3); plot(ah,x1,'sf','byperiod');
title('\rm sf by periods');

ah = subplot(6,3,4); plot(ah,xs,'tr','sa','combined',ticks{:});
ylim([0,12]); ylabel('seas'); title('');
ah = subplot(6,3,5); plot(ah,xs,'sf','ir','combined',ticks{:}); title('');
ah = subplot(6,3,6); plot(ah,xs,'sf','byperiod'); title('');

ah = subplot(6,3,7); plot(ah,x1,'d12','d11','combined',ticks{:});
ylim([0,12]); ylabel('simplified X-11'); title('');
ah = subplot(6,3,8); plot(ah,x1,'d10','d13','combined',ticks{:}); title('');
ah = subplot(6,3,9); plot(ah,x1,'d10','byperiod'); title('');

ah = subplot(6,3,10); plot(ah,x3,'d12','d11','combined',ticks{:});
ylim([0,12]); ylabel('X-13ARIMA'); title('');
ah = subplot(6,3,11); plot(ah,x3,'d10','d13','combined',ticks{:}); title('');
ah = subplot(6,3,12); plot(ah,x3,'d10','byperiod'); title('');

ah = subplot(6,3,13); plot(ah,x3s,'s12','s11','combined',ticks{:});
ylim([0,12]); ylabel('TRAMO-SEATS'); title('');
ah = subplot(6,3,14); plot(ah,x3s,'s10','s13','combined',ticks{:}); title('');
ah = subplot(6,3,15); plot(ah,x3s,'s10','byperiod'); title('');

ah = subplot(6,3,16); plot(ah,x1,'crp','csa', ...
    'combined','selection',[0 0 0 1 0 0 0 0 0 0 0 0 0],ticks{:});
ylim([0,12]); ylabel('camplet'); title('');
ah = subplot(6,3,17); plot(ah,x1,'csf','crp', ...
    'combined','selection',[0 1 0 0 0 0 0 0 0 0 0 0 0],ticks{:});
title('\rm sf and ir'); title('');
ah = subplot(6,3,18); plot(ah,x1,'csf','byperiod');
title('\rm sf by periods'); title('');

%% conclusion

disp(WrapLines(['X-11 is not the original version of the X-11 algorithm. ', ...
    'But the seasonal factors are comparable to X-13 (except for the ', ...
    'edge of the sample, because the simplified X-11 does not use ARIMA ', ...
    'to extend the data). The same is true for seas.'],63));

disp(WrapLines(['A difference, however, is the irregular component. ', ...
    'This is supposed to contain only high frequencies, but in X-11, ', ...
    'seas, and fixedseas, the irregular contains some lower frequency ', ...
    'as well.'],63));

disp(WrapLines(['camplet also produces quite a different result. ', ...
    'First of all, the trend is backward looking and has therefore a ', ...
    'phase shift. Moreover, the seasonal factors are quite volatile with ', ...
    'this algorithm.'],63));
