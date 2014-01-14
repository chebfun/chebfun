function varargout = chebpolyplot(f, varargin)
%CHEBPOLYPLOT   Display Chebyshev coefficients graphically.
%   CHEBPOLYPLOT(F) plots the Chebyshev coefficients of a CHEBFUN F on a
%   semilogy scale. A horizontal line at the epslevel of F is also plotted. If
%   F is an array-valued CHEBFUN or has breakpoints, then a curve is plotted
%   for each FUN of each component (column) of F.
%
%   CHEBPOLYPLOT(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale. If S
%   contains a string 'NOEPSLEVEL' the epslevel is not plotted.
%
%   H = CHEBPOLYPLOT(F) returns a column vector of handles to lineseries
%   objects. The final entry is that of the epslevel plot.
%
% See also CHEBFUN/PLOT

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

% Store the hold state of the current axis:
holdState = ishold;

% Grab some colours:
col = [];
if ( nargin > 1 && ischar(varargin{1}) && numel(varargin{1}) < 4 )
    % (Note. LineStyle definitions have at most 4 characters. Other plotting
    % parameters have more than characters.)
    str = varargin{1};
    col = regexp(varargin{1}, '[bgrcmykw]', 'match');
    col = col{:};
    varargin(1) = [];
end
if ( isempty(col) )
    colIdx = find(cellfun(@(arg) all(ischar(arg)) && any(strfind(arg, 'color')), varargin));
    if ( colIdx )
        col = varargin{colIdx+1};
        varargin(colIdx+(0:1)) = [];
    else
        col = get(0, 'DefaultAxesColorOrder');
    end
end

% Initialise the output:
h1 = cell(numel(f), 1);
h2 = cell(numel(f), 1);

% Loop over the columns:
for k = 1:numel(f)
    % Extract the kth column:
    fk = f(k);
    
    % Colours for this column:
    if ( ischar(col) )
        colk = col;
    else
        numColsFk = numColumns(fk);
        colk = col(1:numColsFk, :);
        col(1:numColsFk, :) = [];
    end

    % Call the column version:
    [h1{k}, h2{k}] = columnChebpolyplot(fk, colk, varargin{:});
    
end

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output if one was requested:
if ( nargout > 0 )
    varargout = {cell2mat(h1), cell2mat(h2)};
end

end

function [h1, h2] = columnChebpolyplot(f, col, varargin)

% Get the coeff data:
numFuns = numel(f.funs);
n = cell(numFuns, 1);
c = get(f, 'coeffs');
% [TODO]: GET('coeffs') should not be used. We proably need chebpolyplotData()
% methods, similar to those used by the PLOT() funtion.
if ( ~iscell(c) )
    c = {c};
end
c = cellfun(@abs, c, 'UniformOutput', false);
for k = 1:numFuns
    n{k} = (size(c{k}, 1)-1:-1:0).';
end
numCols = size(c{1}, 2);

% Get vscale, epslevel data:
v = get(f, 'vscale-local');
e = get(f, 'epslevel-local');
ve = bsxfun(@times, v, e);

% TODO: fix this!
% % Add a tiny amount to zeros to make plots look nicer:
% minve = min(ve, [], 1)
% for k = 1:numfuns
%     % Use smaller of min. of (vscale)*(epslevel) and the smallest nonzero coeff.
%     c{k}(~c{k}) = min(minve(k), min(c{k}(logical(c{k}(:)))));
% end

% Shape it:
data = reshape([n c]', 1, 2*numFuns);

% Plot the coeffs:
h1 = semilogy(data{:}, varargin{:}, 'color', col);
hold on

% Reshape data for epslevel plot:
n = cellfun(@(x) x([end ; 1]), n, 'UniformOutput', false);
n = reshape(repmat(n', numCols, 1), numCols*numel(n), 1);
ve = reshape(ve, numFuns*numCols, 1);
ve = mat2cell(repmat(ve, 1, 2), ones(numFuns*numCols, 1), 2);
data = reshape([n ve]', 1, 2*numFuns*numCols);

% Plot the epslevels:
h2 = semilogy(data{:}, varargin{:}, 'color', col);
for k = 1:numel(h2)
    c = get(h1(k), 'color');
    set(h2(k), 'linestyle', ':', 'linewidth', 1, 'marker', 'none', 'color', c);
end

end

