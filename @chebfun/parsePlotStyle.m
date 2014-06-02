function [lineStyle, pointStyle, jumpStyle, out] = parsePlotStyle(varargin)
%PARSEPLOTSTYLE   Parse inputs to PLOT. Extract 'lineWidth', etc.
%   [L, P, J, OTHER] = PARSEPLOTSTYLE(VARARGIN) parses the inputs VARARGIN and
%   strips out inputs to the MATLAB/PLOT() that should only be in cluded once.
%   For example, 'LineWidth' or 'MarkerSize'. Those options which correspond to
%   the Line part of the plot are returned in L, those corresponding to the
%   discrete points are returned in P, options for the 'jumpLine' appear in J,
%   and all other inputs are returned in OTHER as a cell array.
%
% See also PLOT, PLOT3.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

lineOpts = {'LineStyle', 'LineWidth'};
pointOpts = {'Marker', 'MarkerSize', 'MarkerFaceColor', 'MarkerEdgeColor'};

% Supress MLINT warning for growing arrays.
%#ok<*AGROW>

% Initialise:
lineStyle = {};
pointStyle = {};

% Do JumpLine first.
[jumpStyle, varargin] = parseJumpStyle(varargin{:});

k = 1; % Look at all remaining arguments.
while ( k < numel(varargin) )
    vk = varargin{k};

    if ( any(strcmpi(vk, lineOpts)) )
        % Line option:
        lineStyle = [lineStyle, vk, varargin{k+1}]; 
        varargin(k:k+1) = [];
        
    elseif ( any(strcmpi(vk, pointOpts)) )
        % Point option:
        pointStyle = [pointStyle, vk, varargin{k+1}];
        varargin(k:k+1) = [];
        
    elseif ( strcmpi(vk, 'color') )
        % Option for all:
        lineStyle = [lineStyle, vk, varargin{k+1}];
        pointStyle = [pointStyle, vk, varargin{k+1}];
        jumpStyle = [jumpStyle, vk, varargin{k+1}];
        varargin(k:k+1) = [];
    else
        k = k + 1;
    end

end

% Assign the remaining arguments to OUT:
out = varargin;

end

function [jumpStyle, varargin] = parseJumpStyle(varargin)
%PARSEJUMPSTYLE   Parse the 'jumpline' style for CHEBFUN plot functions.
%   [JUMPSTYLE, VARARGIN] = PARSEJUMPSTYLE(VARARGIN) takes the VARARGIN input
%   for a CHEBFUN plotting command and parses out the 'jumpline' option,
%   converting it into a sequence of name-value pairs suitable for passing to
%   MATLAB's built-in plotting functions, which are stored in the cell array
%   JUMPSTYLE.  The remainder of VARARGIN after removing the 'jumpline' option
%   and its value is returned in the VARARGIN output.

jumpStyle = {};
for k = 1:numel(varargin)
    
    if ( ~strcmpi(varargin{k}, 'jumpline') )
        % We're only looking for the 'jumpLine' flag here.
        continue
    end

    tmp = varargin{k+1};
    varargin(k:(k+1)) = [];
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

    % We only expect one 'jumpLine' input, so we can return now.
    return
end

end
