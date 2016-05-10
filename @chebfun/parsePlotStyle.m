function [lineStyle, pointStyle, jumpStyle, deltaStyle, out] = parsePlotStyle(varargin)
%PARSEPLOTSTYLE   Parse inputs to PLOT. Extract 'lineWidth', etc.
%   [L, P, J, D, OTHER] = PARSEPLOTSTYLE(VARARGIN) parses the inputs VARARGIN and
%   strips out inputs to the MATLAB/PLOT() that should only be in cluded once.
%   For example, 'LineWidth' or 'MarkerSize'. Those options which correspond to
%   the Line part of the plot are returned in L, those corresponding to the
%   discrete points are returned in P, options for the 'jumpLine' appear in J,
%   and all other inputs are returned in OTHER as a cell array.
%
% See also PLOT, PLOT3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

lineOpts = {'LineStyle', 'LineWidth'};
pointOpts = {'Marker', 'MarkerSize', 'MarkerFaceColor', 'MarkerEdgeColor'};

% Suppress MLINT warning for growing arrays.
%#ok<*AGROW>

% Initialise:
lineStyle = {};
pointStyle = {};

% Do JumpLine and DeltaLine first.
[jumpStyle, deltaStyle, varargin] = parseLineStyle(varargin{:});

k = 1; % Look at all remaining arguments.
while ( k < numel(varargin) )
    vk = varargin{k};

    % Using strncmp() here is technically incorrect, but it's easier than
    % writing the code to properly do the name matching, and it's highly
    % unlikely that it will ever actually cause us any problems.
    if ( any(strncmpi(vk, lineOpts, 5)) )
        % Line option:
        lineStyle = [lineStyle, vk, varargin{k+1}]; 
        varargin(k:k+1) = [];
        
    elseif ( any(strncmpi(vk, pointOpts, 6)) )
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

function [jumpStyle, deltaStyle, varargin] = parseLineStyle(varargin)
%PARSEJUMPSTYLE   Parse 'jumpline' and 'deltaLine' style for CHEBFUN plot functions.
%   [JUMPSTYLE, DELTASTYLE, VARARGIN] = PARSEJUMPSTYLE(VARARGIN) takes the 
%   VARARGIN input for a CHEBFUN plotting command and parses out the 
%   'jumpline' and 'deltaline' option, converting it into a sequence of 
%   name-value pairs suitable for passing to MATLAB's built-in plotting 
%   functions, which are stored in the cell array JUMPSTYLE and DELTASTYLE.
%   The remainder of VARARGIN, after removing the 'jumpline' and 'deltaline'
%   options and their values, is returned in the VARARGIN output.

jumpStyle = {};
deltaStyle = {};
N = numel(varargin);
k = 1;
while k <= N 
    % Loop and look for 'jumpline' or 'deltaline':
    if ( strcmpi(varargin{k}, 'jumpline') )
        % Parse the style string:
        jumpStyle = parseStyle(varargin{k+1});
        varargin(k:(k+1)) = [];
        N = N - 2;
    elseif ( strcmpi(varargin{k}, 'deltaline') )
        % Parse the style string:
        deltaStyle = parseStyle(varargin{k+1});
        varargin(k:(k+1)) = [];
        N = N - 2;
    end
    k = k + 1;
end
end

function style = parseStyle(styleString)
%PARSESTYLE   Parse style string for plotting.
%   STYLE = PARSESTYLE(STYLESTRING) takes the input STYLESTRING and 
%   converts into a cell array containing a sequence of name-value pairs 
%   suitable for passing to MATLAB's built-in plotting functions.

style = {};
if ( iscell(styleString) )
    cc = regexp(styleString{1},'[bgrcmykw]', 'match');
    if ( ~isempty(cc) )
        % Forgive " 'jumpline', {'b', ...} " by inserting a 'color'.
        style = ['Color', cc, styleString{2:end}];
    else
        style = styleString;
    end
    return
end

ll = regexp(styleString, '[-:.]+','match');           % style
if ( ~isempty(ll) )
    style = [style, 'LineStyle', ll];
end

cc = regexp(styleString,'[bgrcmykw]', 'match');       % color
if ( ~isempty(cc) )
    style = [style, 'Color', cc];
end

mm = regexp(styleString,'[.ox+*sdv^<>ph]', 'match');  % marker
if ( ~isempty(mm) )
    style = [style, 'Marker', mm];
end

if ( any(strcmpi(styleString, {'none', 'off', ''})) ) % off
    style = {'LineStyle', 'none'};
end

end
