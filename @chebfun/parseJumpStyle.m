function [jumpStyle, varargin] = parseJumpStyle(varargin)
%PARSEJUMPSTYLE   Parse the 'jumpline' style for CHEBFUN plot functions.
%   [JUMPSTYLE, VARARGIN] = PARSEJUMPSTYLE(VARARGIN) takes the VARARGIN input
%   for a CHEBFUN plotting command and parses out the 'jumpline' option,
%   converting it into a sequence of name-value pairs suitable for passing to
%   MATLAB's built-in plotting functions, which are stored in the cell array
%   JUMPSTYLE.  The remainder of VARARGIN after removing the 'jumpline' option
%   and its value is returned in the VARARGIN output.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

jumpStyle = {};
for idx = 1:numel(varargin)
    if ( ~strcmpi(varargin{idx}, 'jumpline') )
        continue
    end

    tmp = varargin{idx+1};
    varargin(idx:(idx+1)) = [];
    if ( iscell(tmp) )
        cc = regexp(tmp{1},'[bgrcmykw]', 'match');
        if ( ~isempty(cc) )
            % Forgive " 'jumpline', {'b', ...} " by inserting a 'color'.
            jumpStyle = ['Color', cc, tmp{2:end}];
        else
            jumpStyle = tmp;
        end
        return
    end

    ll = regexp(tmp, '[-:.]+','match');           % style
    if ( ~isempty(ll) )
        jumpStyle = [jumpStyle, 'LineStyle', ll];
    end

    cc = regexp(tmp,'[bgrcmykw]', 'match');       % color
    if ( ~isempty(cc) )
        jumpStyle = [jumpStyle, 'Color', cc];
    end

    mm = regexp(tmp,'[.ox+*sdv^<>ph]', 'match');  % marker
    if ( ~isempty(mm) )
        jumpStyle = [jumpStyle, 'Marker', mm];
    end

    if ( any(strcmpi(tmp, {'none', 'off', ''})) ) % off
        jumpStyle = {'LineStyle', 'none'};
    end

    return
end

end
