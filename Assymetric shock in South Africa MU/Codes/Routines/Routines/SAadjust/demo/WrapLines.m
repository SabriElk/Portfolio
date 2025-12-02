% Wrap lines at spaces so that they have at most l character per line.
% Preappend leadText to every line.
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% PROBLEM:
% s = ['123',char(10),'123456789',char(10),'123456']
% WrapLines(s,5,'>')
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% eol is supposed to be char(10).
% The platform-specific eol can be gotten with
% double(char(java.lang.System.getProperty('line.separator')))

function str = WrapLines(str,l,leadText)
    if nargin < 3
        leadText = '';
    end
    if nargin < 2
        l = 79;
    end
    str = [leadText,str];
    if ~strcmp(str(end),char(10))
        str = [str,char(10)];
    end
    posLF    = [0,strfind(str,char(10)),length(str)];
    startpos = posLF(find(diff(posLF) > l)); %#ok<*FNDSB>
    while ~isempty(startpos)
        posSP = find(ismember(str(startpos(1)+1:startpos(1)+1+l),' '), ...
            1, 'last') + startpos(1);
        if isempty(posSP)
            % no space available; cut in the middle of a word
            str = [str(1:startpos(1)+l), char(10), ...
                leadText, str(startpos(1)+l+1:end)];
        else
            % replace last available space with lf
            str = [str(1:posSP-1), char(10), ...
                leadText, str(posSP+1:end)];
        end
        posLF    = [1,strfind(str,char(10)),length(str)];
        startpos = posLF(find(diff(posLF) > l+1));
    end
end
