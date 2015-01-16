function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display coefficients graphically.
%   PLOTCOEFFS(F) plots the coefficients underlying the representation of a
%   CHEBFUN F on a semilogy scale. A horizontal line at the epslevel of F is
%   also plotted. If F is an array-valued CHEBFUN or has breakpoints, then a
%   curve is plotted for each FUN of each component (column) of F.
%
%   PLOTCOEFFS(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, as in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale. If S
%   contains a string 'NOEPSLEVEL', the epslevel is not plotted.
%
%   H = PLOTCOEFFS(F) returns a column vector of handles to lineseries objects.
%   The final entry is that of the epslevel plot.
%
%   What 'coefficients' means in this context is dictated by the 'tech' that is
%   being used to represent F. See the PLOTCOEFFS method at the relevant tech
%   level for more details. Note that you can find what tech is being used by
%   calling GET(F, 'tech').
%
% See also CHEBFUN/PLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

% We can only plot the coefficients of one CHEBFUN at a time:
if ( any(cellfun(@(f) isa(f, 'chebfun'), varargin)) )
    error('CHEBFUN:CHEBFUN:plotcoeffs:multipleChebfuns', ...
        'Calls of the form PLOTCOEFFS(F, ''b'',  G, ''r'') are not supported.');
end

% Store the hold state of the current axis:
holdState = ishold;

% Grab some colours:
col = [];
if ( nargin > 1 && ischar(varargin{1}) && numel(varargin{1}) < 4 )
    % (Note. LineStyle definitions have at most 4 characters. Other plotting
    % parameters have more than characters.)
    col = regexp(varargin{1}, '[bgrcmykw]', 'match');
    if ( numel(col) > 1 )
        error('CHEBFUN:plotcoeffs:color', ...
            'Error in color/linetype argument.');
    elseif ( ~isempty(col) )
        col = col{:};
%         varargin(1) = []; % Don't remove as we need marker information.
    end
end
if ( isempty(col) )
    colIdx = find(cellfun(@(arg) all(ischar(arg)) && ...
        any(strfind(arg, 'color')), varargin));
    if ( colIdx )
        col = varargin{colIdx+1};
        varargin(colIdx+(0:1)) = [];
    else
        col = get(0, 'DefaultAxesColorOrder');
        col = repmat(col, ceil(numColumns(f)/7), 1);
    end
end

% Deal with 'LogLog' and 'noEpsLevel' input:
doLogLog = cellfun(@(s) strcmpi(s, 'loglog'), varargin);
% varargin(doLoglog) = [];
doLogLog = any(doLogLog);
noEpsLevel = cellfun(@(s) strcmpi(s, 'noEpsLevel'), varargin);
varargin(noEpsLevel) = []; % Strip this out, as we don't want to pass it down.
doEpsLevel = ~any(noEpsLevel);

% Convert to a cell array for easy handling:
f = mat2cell(f);

% Initialise the output:
h1 = cell(1, numel(f));
h2 = cell(1, numel(f));

% Loop over the columns in a quasimatrix:
for k = 1:numel(f)
    % Colours for this columns:
    if ( ischar(col) )
        colk = col;
    else
        colk = col(k,:);
    end
    % Call the column version:
    [h1{k}, h2{k}] = columnPlotCoeffs(f{k}, colk, varargin{:});
end

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Set xScale to logarithmic if requested:
if ( doLogLog )
    set(gca, 'xScale', 'Log');
end

% Remove epslevels if they're not wanted:
if ( ~doEpsLevel )
    for k = 1:numel(f)
        set(h2{k}, 'visible', 'off');
    end
end

% Give an output if one was requested:
if ( nargout > 0 )
    try 
        h1 = cell2mat(h1);
        h2 = cell2mat(h2);
    catch
        % shrug
    end
    varargout = {h1, h2};
end

end

function [h1, h2] = columnPlotCoeffs(f, col, varargin)
numFuns = numel(f.funs);

% Initialise handle storage:
h1 = zeros(numFuns, 1);
h2 = zeros(numFuns, 1);

% Call plotcoeffs at the tech level:
for j = 1:numFuns
    [h1(j), h2(j)] = plotcoeffs(f.funs{j}, varargin{:}, 'color', col); hold on
end

end

