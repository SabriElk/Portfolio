% X11 computes a simplified X-11 seasonal adjustment. It is a bit less refined
% than the Census Bureau original, but it has the advantage that it works with
% arbitrary frequencies, not just monthly or quarterly.
%
% CAUTION: The program computes only a simplified version of the original X-11
% algorithm. Most notably, there is no account of trading day corrections, and
% extreme values are also not taken care of. Still, the adjustment is often
% quite good, the program is applicable to any seasonality (not just 4 and 12 as
% the US Census Bureau program), and the code of the program can be instructive
% to develop variations.
% Later versions of X-11 used an estimated ARIMA model to produce fore- and
% backcasts in order to alleviate the edge of sample problem that offurs in any
% filtering using moving averages. This program does not use ARIMA, but instead
% 'mirrors' at the left and right and applies the filtering after that. This
% simple technique appears to get rid of the edge of sample problem rather well.
%
% NOTE: This program does *not* use the original X-11 executable program from
% the US Census Bureau. It does not support many of the options of that program
% either. This program merely tries to replicate the most important steps of the
% seasonal adjustment performed by the X-11 algorithm using Matlab directly. In
% other words, this is a Matlab implementation of a somewhat simplified version
% of the X-11 algorithm.
% This fact also implies that, unlike the U.S: Census software, this
% implementation accommodates arbitrary frequencies, not just monthly or
% quarterly. This program is just a small addition to the toolbox that makes it
% more complete.
%
% Usage:
%   s = x11(data,period);
%   s = x11([dates,data],period);
%   s = x11(... ,transform);
%
%   s   This is a structure containing the following components:
%       s.p     = period
%       s.type  = type of decomposition
%       s.dates = dates vector
%       s.dat   = data vector
%       s.d13   = trend
%       s.d11   = seasonally adjusted data
%       s.d10   = seasonal factor (cycle)
%       s.d13   = irregular component
%               .
%               .
%               .
%       The other components are from intermediate computation steps. Their
%       meaning can be taken from the documentation of x13as.exe.
%
% transform  Must be one of the following: 'additive','none','multiplicative',
%       'logadditive','pseudoadditive'. It indicates the type of decomposition.
%       'additive' or 'none' : data = tr + sf + ir, sa = tr + ir.
%           'multiplocative' : data = tr * sf * ir, sa = tr * ir.
%              'logadditive' : log(data) = log(tr*sf*ir), sa = exp(tr+ir).
%           'pseudoadditive' : data = tr * (sf + ir - 1), sa = tr * ir.
%
% REMARK: This program uses several smaller programs (trendfilter, seasfilter,
% normalize) that can be used to create a custom seasonal adjustment algorithm
% relatively easily. To understand how, just study the source code of this
% program.
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
% 2017-09-19    Version 1.33    First version of X-11 implementation.

function s = x11(data,p,varargin)

    % validate parameters
    assert(isnumeric(p) && fix(p) == p && numel(p) == 1 && p > 0, ...
        'X13TBX:x11:IllegalPeriod', ...
        'The second argument of x11 must be a positive integer.');
    [row,col] = size(data);
    assert(row>1 && (col == 1 || col == 2), ...
        'X13TBX:x11:NoVector', 'x11 expects a vector, but you have provided a %ix%i array.', ...
            row, col);
    
    % parse args
    if col > 1
        dates = data(:,1);
        data  = data(:,2:end);
    else
        dates = (1:size(data,1))';
    end
    
    % default
    ismult   = true;
    islogadd = false;
    ispseudo = false;
    type     = 'multiplicative';
    
    % interpret varargin
    legal = {'additive','none','multiplicative', ...
        'logadditive','pseudoadditive'};
    
    while ~isempty(varargin)
        type = validatestring(varargin{1},legal);
        switch type
            case {'additive','none'}
                ismult   = false;
                islogadd = false;
                ispseudo = false;
                type = 'additive';
            case 'multiplicative'
                ismult   = true;
                islogadd = false;
                ispseudo = false;
            case 'logadditive'
                ismult   = false;
                islogadd = true;
                ispseudo = false;
            case 'pseudoadditive'
                ismult   = true;
                islogadd = false;
                ispseudo = true;
        end
        varargin(1) = [];
    end
    
    % log-additive ?
    if islogadd
        assert(all(data(~isnan(data)) > 0), ...
            'X13TBX:x11:NegLog', ['Data must be strictly positive ', ...
            'for log-additive decomposition.']);
        data = log(data);
    end
    
    % do the twist
    
    % - stage 1
    tr_1 = trendfilter(data,'cma',p,'mirror',p);
    si_1 = normalize(data,tr_1,ismult);
    sf_1 = seasfilter(si_1,p,'ma',[3,3]);
    if ispseudo
        ir_1 = normalize(si_1,sf_1,false) + 1;
        sa_1 = tr_1 .* ir_1;
    else
        ir_1 = normalize(si_1,sf_1,ismult);
        sa_1 = normalize(data,sf_1,ismult);
    end
    
    % - stage 2
    tr_2 = trendfilter(sa_1,'cma',p,'mirror',p);
    si_2 = normalize(data,tr_2,ismult);
    sf_2 = seasfilter(si_2,p,'ma',[3,3]);
    if ispseudo
        ir_2 = normalize(si_2,sf_2,false) + 1;
        sa_2 = tr_2 .* ir_2;
    else
        ir_2 = normalize(si_2,sf_2,ismult);
        sa_2 = normalize(data,sf_2,ismult);
    end
    
    % - stage 3
    tr_3 = trendfilter(sa_2,'cma',p,'mirror',p);
    si   = normalize(data,tr_3,ismult);
    sf   = seasfilter(si,p,'ma',[5,5]);
    if ispseudo
        ir = normalize(si,sf,false) + 1;
        sa_3 = tr_3 .* ir;
    else
        ir = normalize(si,sf,ismult);
        sa_3 = normalize(data,sf,ismult);
    end
    tr = trendfilter(sa_3,'spencer','mirror',p);
    if ismult
        sa = tr .* ir;
    else
        sa = tr + ir;
    end
    
    % un-log ?
    if islogadd
        data = exp(data);
        tr_1 = exp(tr_1);
        tr_2 = exp(tr_2);
        tr_3 = exp(tr_3);
        tr   = exp(tr);
        si_1 = exp(si_1);
        si_2 = exp(si_2);
        si   = exp(si);
        sf_1 = exp(sf_1);
        sf_2 = exp(sf_2);
        sf   = exp(sf);
        ir_1 = exp(ir_1);
        ir_2 = exp(ir_2);
        ir   = exp(ir);
        sa_1 = exp(sa_1);
        sa_2 = exp(sa_2);
        sa_3 = exp(sa_3);
        sa   = exp(sa);
    end    
    
    % collect everything
    s = struct(...
        'period',   p,      ...
        'type',     type,   ...
        'dates',    dates,	...
        'dat',      data,   ...
        'd12',      tr,     ...
        'd11',      sa,     ...
        'd10',      sf,     ...
        'd8',       si,     ...
        'd13',      ir,     ...
        'd7',       tr_3,   ...
        'd6',       sa_3,   ...
        'c10',      sf_2,   ...
        'c11',      sa_2,   ...
        'c13',      ir_2,   ...
        'c7' ,      tr_2,   ...
        'c9' ,      si_2,   ...
        'b10',      sf_1, 	...
        'b11',      sa_1, 	...
        'b13',      ir_1, 	...
        'b7' ,      tr_1, 	...
        'b8' ,      si_1);
    
end
