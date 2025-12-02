% FIXEDSEAS computes a simple seasonal filter with fixed seasonal factors.
%
% Usage:
%   s = fixedseas(data,period);
%   s = fixedseas([dates,data],period);
%   s = fixedseas(... ,transform);
%   s = fixedseas(... ,method);
%   s = fixedseas(... ,method,methodarg);
%   [s,aggr] = fixedseas(...);
%
% data must be a vector. fixedseas is NaN tolerant, meaning data can contain
% NaNs.
%
% period is a positive number which indicates the length of the seasonal
% cycle (i.e. period = 12 for monthly data, period = 7 for daily data
% having a weekly cycle, or period = 5 if the data is weekdaily).
% 
% The optional arguments determine if the filtering should be done
% additively or multiplicatively, and the method of filter to use for
% computing the trend.
%
% 'transform' is one of the following:
%   'none' or 'additive'        The decomposition is done additively. This
%                               is the default.
%   'logadditive'               The log is applied to the data, the
%                               decomposition is then applied additively, and
%                               the exponential of the result is returned.
%   'multiplicative'            The decomposition is done multiplicatively.
%
% 'method' (and 'methodarg') determines the method of trend. For possible
% choices, see help trendfilter. Default is a centered moving average with
% length equal to period ('cma',period).
%
% 'period' can also be a positive vector. In that case, the seasonal filtering
% is performed several times, removing cycles at all desired frequencies. In
% that case, 'transfrom' and 'method' can be cellarrays, containing one method
% (plus argument) for each period. The returned s is then a structure with as
% many components as there are components in 'period'.
%
% If period is a vector and transfrom is the same for each period, an additional
% aggregated structure is appended to s, providing the cumulated seasonal
% factors etc. In that case, aggr is returned as true.
%
% The program will select default values for 'lambda','roughness', or
% 'degree', respectively, if you do not specify them. If you use a vector
% for the 'period' argument (filtering out multiple periods), then you can
% also specify vectors of lambda/roughness/degree-arguments, one for each
% component of your period-vector.
%
% s is a struct with the following fields:
%   .period     Period(s) that has/have been filtered.
%   .transform  Either 'none' or 'log' or 'mult'.
%   .method     The method used for computing the trend.
%   .methodarg  possibly a parameter for the smoothing algorithm.
%   .dates      The original dates. If none were provided, this is just a vecor
%               counting from 1 to the number of data points.
%   .dat        The original data.
%   .tr         Long term trend (by default the moving average, but other
%               choices are possible, see above).
%   .sa         Seasonally adjusted series (= dat – sf, or exp(dat – sf),
%               respectively).
%   .sf         Seasonal factors.
%   .ir         Irregular (= sa – tr or exp(sa – tr), respectively).
%
% Data is decomposed into the three components, trend (tr), seasonal factor
% (sf), and irregular (ir). For the additive decomposition, it is always
% the case that data = tr + sf + ir. Furthermore, sa = data - sf (or
% equivalently, sa = tr + ir). For the multiplicative decomposition, data =
% tr * sf * ir, and sa = data ./ sf (or equivalently, sa = tr * ir).
%
% Example 1:
%   truetrend = 0.02*(1:200)' + 5;
%   % truecycle = sin((1:200)'*(2*pi)/20);
%   truecycle = repmat([zeros(7,1);-0.6;zeros(11,1);0.9],ceil(200/20),1);
%   truecycle = truecycle(1:200);
%   truecycle = truecycle - mean(truecycle);
%   trueresid = 0.2*randn(200,1);
%   data = truetrend + truecycle + trueresid;
%   s = fixedseas(data,20);
%   figure('Position',[78 183 505 679]);
%   subplot(3,1,1); plot([s.dat,s.sa,s.tr,truetrend]); grid on;
%   title('unadjusted and seasonally adjusted data, estimated and true trend')
%   subplot(3,1,2); plot([s.sf,truecycle]); grid on;
%   title('estimated and true seasonal factor')
%   subplot(3,1,3); plot([s.ir,trueresid]); grid on;
%   title('estimated and true irregular')
%   legend('estimated','true values');
%
% Example 2 (multiple cycles):
%   truecycle2 = 0.7 * sin((1:200)'*(2*pi)/14);
%   data = truetrend + truecycle + truecycle2 + trueresid;
%   s = fixedseas(data,[14,20],'hp');
%   figure('Position',[78 183 505 679]);
%   subplot(3,1,1); plot([s.dat,s.sa,s.tr,truetrend]); grid on;
%   title('unadjusted and seasonally adjusted data, estimated and true trend')
%   subplot(3,1,2); plot([s.sf,truecycle+truecycle2]); grid on;
%   title('estimated and true seasonal factor')
%   subplot(3,1,3); plot([s.ir,trueresid]); grid on;
%   title('estimated and true irregular')
%   legend('estimated','true values');
%
% Note that fixedseas(data,[14,20]) is not the same as
% fixedseas(data,[20,14]). The filters are applied iteratively, from left
% to right. The ordering matters, so the results differ.
%
% Detailed description of the model: Let x be some timeseries. As an
% example, we compute fixedseas(x,6).
% *** STEP 1 ***
% We compute a 6-period centered moving average,
%   trend(t) = sum(0.5x(t-3)+x(t-2)+x(t-1)+x(t)+x(t+1)+x(t+2)+0.5x(t+3))/6
% The weights on the extreme values of the window are adapted so that the
% sum of the weights is equal to period. So, for instance, if period = 7,
% the weight on x(t-3) and x(t+3) would be 1.0; if period = 6.5, the weight
% would be 0.75.
% [Note: By default the trend is computed as the centered moving average,
% and this is what is explained here. Other specifications are possible,
% namely detrend, hodrick-prescott, spline, and polynomial.]
% *** STEP 2 ***
% Compute the individual deviations of x from the trend,
%   d = x - trend.
% *** STEP 3 ***
% Compute the average deviation over all observations on a cycle of 6
% periods,
%   m(1) = mean(d(1) + d(7) + d(13) + d(19) + ...)
%   m(2) = mean(d(2) + d(8) + d(14) + d(20) + ...)
%   ...
%   m(6) = mean(d(6) + d(12) + d(18) + d(24) + ...)
% *** STEP 4 ***
% Normalize m so that its average is zero,
%   n = (m(1)+m(2)+...+m(6))/6
%   sf(1) = m(1) - n, sf(2) = m(2) - n, ..., sf(6) = m(6) - n
% These are the seasonal factors.
% *** STEP 5 ***
% Compute the seasonally adjusted time series as sa = x - sf.
% *** STEP 6 ***
% Compute the irregular as ir = sa - trend. This is the part of the
% fluctuations of x that is not explained by the seasonal factors or the
% trend (= moving average).
%
% STEP 1 as described here is for the 'moving averqage' trend type, which
% is the default. This step is different for the different trend types
% that are available. STEP 2 to 6 are, however, independent of the type of
% trend that is computed.
%
% If the multiplicative option is used, the logarithm of the data is
% processed and the exponential of the processed time series is returned.
% So, s = fixedseas(data,period,'log') is materially the same as
% s2 = fixedseas(log(data),period). Then, exp(s2.sa) = s.sa,
% exp(s2.sf) = s.sf, and exp(s2.tr) = s.tr.
%
% NOTE: This program is part of the X-13 toolbox, but it is completely
% independent of the Census X-13 program. It uses a much simpler strategy
% to filter seasonal cycles, but it has the advantage of accommodating arbitrary
% frequencies, not just monthly or quarterly as X-13 does (the SEATS
% specification supports a few more frequencies). This program is just a small
% addition to the toolbox that makes it more complete.
%
% see also guix, x13, makespec, x13spec, x13series, x13composite, 
% x13series.plot,x13composite.plot, x13series.seasbreaks,
% x13composite.seasbreaks, fixedseas, camplet, spr, InstallMissingCensusProgram
%
% Author  : Yvan Lengwiler
% Version : 1.33
%
% If you use this software for your publications, please reference it as:
%
% Yvan Lengwiler, 'X-13 Toolbox for Matlab, Version 1.33', Mathworks File
% Exchange, 2014-2018.
% url: https://ch.mathworks.com/matlabcentral/fileexchange/49120-x-13-toolbox-for-seasonal-filtering

% History:
% 2018-09-19    Version 1.33    New implementation of fixedseas using the
%                               programs in seas.
% 2016-09-06    Version 1.18.4  breakpoints of detrend method are now given as
%                               datevectors (not as indexes indicating the
%                               positions of the breaks) if dates are given by
%                               the user.
% 2016-09-05    Version 1.18.3  Multiple tr, ir, and si when multiple periods
%                               are selected.
% 2016-07-10    Version 1.17.1  Improved guix. Bug fix in x13series relating to
%                               fixedseas.
% 2016-07-06    Version 1.17    First release featuring guix. Bug fix in the
%                               computation of 'ir' when decomposition is
%                               multiplicative.
% 2016-03-03    Version 1.16    Adapted to X-13 Version 1.1 Build 26.
% 2015-08-20    Version 1.15    Significant speed improvement. The imported
%                               time series will now be mapped to the first
%                               day of month if this is the case for the
%                               original data as well. Otherwise, they will
%                               be mapped to the last day of the month. Two
%                               new options --- 'spline' and 'polynomial'
%                               --- for fixedseas. Improvement of .arima,
%                               bugfix in .isLog.
% 2015-08-14    Version 1.14.2  Added 'spline' and 'polynomial' trend
%                               types. Added default typearg values for all
%                               trend types.
% 2015-07-25    Version 1.14    Improved backward compatibility. Overloaded
%                               version of seasbreaks for x13composite. New
%                               x13series.isLog property. Several smaller
%                               bugfixes and improvements.
% 2015-07-24    Version 1.13.3  Resolved some backward compatibility
%                               issues (thank you, Carlos).
% 2015-07-07    Version 1.13    seasma removed, replaced by fixedseas.
%                               Complete integration of fixedseas into
%                               x13spec, with fore-/backcast extension
%                               before computing trend for simple seasonal
%                               adjustment. Various improvemnts to
%                               x13series.plot (including 'separate' 
%                               option). seasbreaks program to identify
%                               seasonal breaks. Better support for PICKMDL
%                               model list files. Added '-n' to list of
%                               default flags in x13. Select print requests
%                               added as default in makespec.
% 2015-05-21    Version 1.12    Several improvements: Ensuring backward
%                               compatibility back to 2012b (possibly
%                               farther); Added 'fixedseas' option to x13;
%                               Added Runsfixedseas to x13series; other
%                               improvements throughout. Changed numbering
%                               of versions to be in synch with FEX's
%                               numbering.
% 2015-05-18    Version 1.6.1   removed epanechnikov option (it was stupid
%                               to begin with)
% 2015-04-28    Version 1.6     x13as V 1.1 B 19
% 2015-03-13    Version 1.4     'detrend' with break points added
% 2015-02-14    Version 1.3     'detrend' and 'HP' are now NaN-tolerant
% 2015-02-03    Version 1.2     support for Epanechnikov, Hodrick-Prescott,
%                               and detrend
% 2015-01-30    Version 1.1     support for fractional period argument and
%                               for multiplicative decomposition
% 2015-01-26    Version 1.0d    residuals called .rsd now
% 2015-01-25    Version 1.0c    improved help
% 2015-01-24    Version 1.0b    small bugfix
% 2015-01-22    Version 1.0     first version

%#ok<*TRYNC>
%#ok<*AGROW>

function [s,aggregate] = fixedseas(data,p,varargin)

    % validate arguments
    if iscell(p); p = cellfun(@(c) str2double(c), p); end
    if ischar(p); p = str2num(p) ; end %#ok<ST2NM>
    assert(isnumeric(p) && all(fix(p) == p) && all(p > 0), ...
        'X13TBX:x11:IllegalPeriod', ['The second argument of x11 must ', ...
        'contain only positive integers.']);
    [rows,cols] = size(data);
    assert(rows>1 && (cols == 1 || cols == 2), ...
        'X13TBX:x11:NoVector', ['x11 expects a vector, but you have ', ...
        'provided a %ix%i array.'], rows, cols);
    
    % parse args
    if size(data,2) > 1
        dates = data(:,1);
        data  = data(:,2:end);
    else
        dates = (1:size(data,1))';
    end
    
    transform = [];
    method    = [];
    methodarg = [];

    kerneltypes = {'rehomme-ladiray','bongard', 'henderson', ...
        'spencer','spencer15','cma','centered moving average', 'ma', ...
        'moving average','uniform','rectangular','rectangle','box', ...
        'epanechnikov','triangle','triangular','biweight','quartic', ...
        'triweight','tricube','cosine','optcosine','gaussian','normal', ...
        'cauchy','picard'};
    legalmethods = {'mean','deviation','reldeviation', 'detrend','spline',...
        'polynomial','hp',kerneltypes{:}};
    legaltransform = {'additive','none','multiplicative','logadditive'};

    possibleMethodarg = false;
    while ~isempty(varargin)
        if iscell(varargin{1})  % a cell in a cell
            arg = varargin{1}{1};
        else
            arg = varargin{1};
        end
        try
            validatestring(arg,legaltransform);
            possibleMethodarg = false;
            if ~isempty(transform)
                warning(['overwriting former transform ''',transform{1},'''']);
            end
            transform{end+1} = varargin{1};
            if ~iscell(transform); transform = {transform}; end
        catch
            try
                validatestring(arg,legalmethods);
                possibleMethodarg = true;
            method{end+1} = varargin{1};
            if ~iscell(method); method = {method}; end
            catch
                if possibleMethodarg
                    methodarg{end+1} = varargin{1};
                else
                    disp(varargin{1});
                    err = MException('X13TBX:fixedseas:ill_arg', ...
                        'There is an illegal argument here.');
                    throw(err);
                end
            end
        end
        varargin(1) = [];
    end
    
    % defaults
    if isempty(transform); transform = {'multiplicative'}; end
    if isempty(method);    method    = {'cma'};            end
    if isempty(methodarg); methodarg = {[]};               end
    
    % get correct number of args; repeat if necessary
    p = num2cell(p); np = numel(p); data = [{data},cell(1,np-1)];
    transform = reparg(transform,np,false);
    method    = reparg(method   ,np,false);
    methodarg = reparg(methodarg,np,true);
    
    % get non-abbreviated versions of parameters
    for q = 1:np
        transform{q} = validatestring(transform{q},legaltransform);
        if strcmp(transform{q},'none'); transform{q} = 'additive'; end
        method{q}    = validatestring(method{q}   ,legalmethods  );
    end

    % do the work (one period at a time)
    
    [tr{1},si{1},sf{1},sa{1},ir{1}] = ...
        seas1p(data{1},p{1},transform{1},method{1},methodarg{1});
    for q = 2:np
        data{q} = sa{q-1};
        [tr{q},si{q},sf{q},sa{q},ir{q}] = ...
            seas1p(data{q},p{q},transform{q},method{q},methodarg{q});
    end
    
    % add a cumulative adjustment if all periods use the same transformation
    
    if np>1 && all(strcmp(transform{1},transform))
        aggregate = true;
        p{end+1}            = p;
        transform{end+1}    = transform{1};
        method{end+1}       = method;
        methodarg{end+1}    = methodarg;
        data{end+1}         = data{1};
        tr{end+1}           = tr{end};
        sa{end+1}           = sa{end};
        if all(strcmp(transform,'multiplicative'))
            si{end+1} = si{1}; sf{end+1} = sf{1}; ir{end+1} = ir{1};
            for r = 2:np
                si{end} = si{end}.* si{r};
                sf{end} = sf{end}.* sf{r};
                ir{end} = ir{end}.* ir{r};
            end
        else
            si{end+1} = si{1}; sf{end+1} = sf{1}; ir{end+1} = ir{1};
            for r = 2:np
                si{end} = si{end} + si{r};
                sf{end} = sf{end} + sf{r};
                ir{end} = ir{end} + ir{r};
            end
        end
    else
        aggregate = false;
    end
    
    % pack everything up
    
    s = struct( ...
        'period',       p,          ...
        'transform',    transform,  ...
        'method',       method,     ...
        'methodarg',    methodarg,  ...
        'dates',        dates,      ...
        'dat',          data,       ...
        'tr',           tr,         ...
        'sa',           sa,         ...
        'sf',           sf,         ...
        'si',           si,         ...
        'ir',           ir);

    % --------------------------------------------------------------------------
    
    % filter one frequency
    function [trNew,siNew,sfNew,saNew,irNew] = seas1p(thisdata,thisp, ...
            thistrans,thismethod,thismethodarg)
        % interpret transform
        switch thistrans
            case 'additive'
                ismult   = false;
                islogadd = false;
            case 'multiplicative'
                ismult   = true;
                islogadd = false;
            case 'logadditive'
                ismult   = false;
                islogadd = true;
                thisdata = log(thisdata);   % take log of data in this case
        end
        % default methodarg if none was given
        if isempty(thismethodarg)
            switch thismethod
                case 'spline'
                    h = ((dates(end)-dates(1)) / (numel(dates)-1)) ./ thisp;
                    thismethodarg = 1 ./ (1 + h.^3 / 0.6);
                case 'polynomial'
                    thismethodarg = floor(numel(dates) ./ thisp);
                case 'hp'
                    slope = 5.91863781313348;
                    absolute = -7.10636;
                    thismethodarg = exp(absolute + slope * log(thisp));
%                     thismethodarg = 8.25 * thisp^3.9;   % see Ravn-Uhlig (http://discovery.ucl.ac.uk/18641/1/18641.pdf)
                case {'rehomme-ladiray','henderson','bongard'}
                    thismethodarg = 2*thisp-1;
                case 'detrend'
                    thismethodarg = [];
                otherwise
                    thismethodarg = thisp;
            end
        end
        % do the computations
        % -- method if detrend (that requires some additional work)
        if strcmp(thismethod,'detrend')
            if ~isempty(thismethodarg)  % detrend with breakpoints
                switch class(thismethodarg)
                    case 'char'
                        thismethodarg = str2num(thismethodarg); %#ok<ST2NM>
                    case 'cell'
                        thismethodarg = cell2mat(thismethodarg);
                end
                breakpoints = NaN(1,numel(thismethodarg));
                for b = 1:numel(thismethodarg)
                    breakpoints(b) = find(dates<=thismethodarg(b),1,'last');
                end
                breakpoints(isempty(breakpoints)|isnan(breakpoints)) = [];
                trNew = trendfilter(thisdata,'detrend',breakpoints, ...
                    'mirror',thisp);
            else
                trNew = trendfilter(thisdata,'detrend', ...
                    'mirror',thisp);
            end
        else
        % -- method is not 'detrend'
            trNew = trendfilter(thisdata,thismethod,thismethodarg, ...
                'mirror',thisp);
        end
        siNew = normalize(thisdata,trNew,ismult);
        if ismult
            sfNew = seasfilter(siNew,thisp,'reldeviation');
        else
            sfNew = seasfilter(siNew,thisp,'deviation');
        end
        saNew = normalize(thisdata,sfNew,ismult);
        irNew = normalize(saNew,trNew,ismult);
        % exponentiate if log was taken before
        if islogadd
            trNew = exp(trNew);
            saNew = exp(saNew);
            sfNew = exp(sfNew);
            siNew = exp(siNew);
            irNew = exp(irNew);
        end
    end

    % repeat argument np times and make a cell array
    function r = reparg(arg,np,evaluate)
        if isnumeric(arg)
            arg = mat2str(arg);
        end
        if evaluate && ischar(arg)
            arg = strrep(arg,'(','{'); arg = strrep(arg,')','}');
            arg = strrep(arg,'[','{'); arg = strrep(arg,']','}');
            arg = strrep(arg,' ',','); arg = strrep(arg,';',',');
            arg = eval(arg);
        end
        if ~iscell(arg); arg = {arg}; end
        if np == 1
            r = arg;
        else
            r = cell(1,np);
            r(1:numel(arg)) = arg;
            r(numel(arg)+1:end) = arg(end);
        end
    end

end
