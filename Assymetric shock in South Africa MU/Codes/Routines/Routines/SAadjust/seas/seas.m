% SEAS is a simple program that demonstrates the use of the programs in the
% seas directory of the X-13 toolbox. These programs are:
%
% trendfilter.m         Computes a trend (i.e. smoothed) version of the data.
%                       You can choose from many different methods to do that.
% seasfilter.m          Splits data using splitperiods, smoothes them with
%                       trendfilter, and joins them together again with
%                       joinperiods.
% normalize.m           Computes the difference or the ratio, respectively, of
%                       two series, depending on whether the decomposition is
%                       additive or multiplicative, respectively.
% splitperiods.m        Splits the data into their periods. For instance, with
%                       monthly data, split periods makes twelve times series
%                       out of your data, one for each month of the year.
% joinperiods.m         Reverse of splitperiods.
% fillholes.m           Linear interpolation of missing values.
% wmean.m               Computes the weighted mean. Similar to Matlab's conv
%                       command, but with smarter treatment of the edge of the
%                       data.
% kernelweights.m       Computes the weights for a wide range of kernels. Used
%                       by trendfilter.m and seasfilter.m in conjunction with
%                       wmean.m
% fixedseas.m           A rather elaborate that produces a rather simple version
%                       of seasonal adjustment in which the seasonal factors are
%                       kept fixed over the years.
% x11.m                 An simplementation of a much simplified version of the
%                       original X-11 method of the U.S. Census Bureau.
% camplet.m             A form of seasonal adjustment that was recently
%                       developed and that does not produce revisions when data
%                       are added to the time series. It does that because the
%                       smoothing is completely backward looking (so no centered
%                       filters at all). camplet is separately implemented and
%                       does not use the other tools provided here.
% seas.m                This file. You can experiment with the implementation
%                       in this file, and develop your own seasonal adjustment
%                       routine starting from seas.m.
%
% Usage of seas.m
%   s = seas(data,p)
%   s = seas(data,p,type)
%   [s,x] = seas(...)
%
% data is either a column vector or an array with two columns. In that case,
% the left column is a date vactor and the right column is the data vector.
%
% p is the period of the data that is to be filtered out. So, typically, with
% monthly data, for instance, p should be set to 12.
%
% type is either 'additive' or 'multiplicative'. 'additive' implies that the
% seasonal factor will be zero on average and are subtracted from the unadjusted
% data to get to the seasonally adjusted data. If type is 'multiplicative' the
% seasonal factor is one on average, and the data is divided by the seasonal
% factor to get the seasonally adjusted data.
%
% s contains the output in a struct, x contains that information but packed
% into an x13series object.
%
% NOTE: This program is part of the X-13 toolbox, but it is completely
% independent of the Census X-13 program. It is part of the 'seas' addition to
% the toolbox which allows to implement seasonal filters without using the
% Census Bureau programs.
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
% 2018-09-19    Version 1.33    First version of the 'seas' part of the X-13
%                               toolbox.

function [s,x] = seas(data,p,type)

    % PARSE ARGS
    
    % -- dates present?
    nobs = size(data,1);
    if size(data,2) > 1
        dates = data(:,1);
        data  = data(:,2:end);
    else
        dates = (1:nobs)';
    end
    
    % -- additive or multiplicative?
    if nargin<3 || isempty(type) || all(isnan(type))
        type = 'multiplicative';
    end
    
    type = validatestring(type,{'additive','multiplicative'});
    ismult = strcmp(type,'multiplicative');
    
    % DO THE WORK
    % There are essentially two places where you need to make a choice, and
    % these choices will determine the shape of the seasonal decomposition.
    % The first is: How do you want to smooth the data to compute the trend.
    % Allowing for a lot  of roughness in the trend produces smaller seasonal
    % factors over all. Enforcing a very smooth trend, on the other hand, will
    % produce a more important role for the seasonal factor.
    % The second decision is how you want to smooth the seasonal deviations
    % over consecutive years. If you allow for a lot of roughness here, the
    % seasonal factors will change quickly from year to year, so that the
    % seasonality might appear rather unstables. On the other hand, you can
    % also keep the seasonal factors completely fixed from one year to the
    % next, but that will come at the cost of increasing the 'unexplained'
    % irregular component and making it more serially correlated.
    
    % -- TR (choice #1)
    tr = trendfilter(data,'cma',p,p/2,'mirror',p);
    % b = exp(-7.10636 + 5.91863781313348 * log(p));
    % tr = trendfilter(data,'hp',b,'mirror',p);
    % tr = trendfilter(data,'spencer','mirror',p);
    % tr = trendfilter(data,'henderson',2*p-1,'mirror',p);
    % tr = trendfilter(data,'bongard',2*p-1,'mirror',p);
    % tr = trendfilter(data,'rehomme-ladiray',[2*p-1,3,0.5],'mirror',p);
    % tr = trendfilter(data,'triangle',2*p-1,'mirror',p);
    % tr = trendfilter(data,'epanech',2*p-1,'mirror',p);
    % h = ((dates(end)-dates(1)) / (numel(dates)-1)) ./ p;
    % roughness = 1 ./ (1 + h.^3 / 0.6);
    % tr = trendfilter(data,'spline',roughness,'mirror',p);
    
    % -- SI
    si = normalize(data,tr,ismult);
    
    % -- SF (choice #2)
    % sf = seasfilter(si,p,'cma',5,4,4,'mirror',ceil(nobs/p));
    % sf = seasfilter(si,p,'hp',1400,'mirror',ceil(nobs/p));
    % sf = seasfilter(si,p,'poly',4,'mirror',ceil(nobs/p));
    sf = seasfilter(si,p,'spline',0.1,'mirror',ceil(nobs/p));
    % sf = seasfilter(si,p,'henderson',13,'mirror',ceil(nobs/p));
    % sf = seasfilter(si,p,'epanech',5,'mirror',ceil(nobs/p));
    %
    % To get seasonal factors that are fixed over the years, you need to smooth
    % this using a simple average over the whole sample. How you take the
    % sample dependy on whether you decompose additively or multiplicatively.
    % if ismult
    %     sf = seasfilter(si,p,'reldeviation');
    % else
    %     sf = seasfilter(si,p,'deviation');
    % end
    
    % -- SA and IR
    sa = normalize(data,sf,ismult);
    ir = normalize(sa,tr,ismult);
    
    % COLLECT EVERYTHING
    
    s = struct(...
        'p',        p,		...
        'type',     type,	...
        'dates',    dates,	...
        'dat',      data,	...
        'tr',       tr,		...
        'sa',       sa,		...
        'sf',       sf,		...
        'ir',       ir,		...
        'si',       si);
    
    % pack everything into an x13series object if the user has requested it
    
    if nargout > 1
    
        x = x13series;
        x.prog = 'seas.m';
        x.progtype = 'Custom seasonal filter';
        x.addvariable('dat',dates,data,'dat',1);
        x.addvariable('tr' ,dates,tr  ,'tr' ,1);
        x.addvariable('sa' ,dates,sa  ,'sa' ,1);
        x.addvariable('sf' ,dates,sf  ,'sf' ,1);
        x.addvariable('si' ,dates,si  ,'si' ,1);
        x.addvariable('ir' ,dates,ir  ,'ir' ,1);
        
        warning_state = warning('off','X13TBX:miss_toolbox');
        x.addspectrum('sa' ,1,'sjs', ...
            'Spectrum of first-differenced adjusted (seas)');
        x.addspectrum('ir' ,0,'sgs', ...
            'Spectrum of irregular (seas)');
        warning(warning_state);
    
    end

end
