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

% Get the coeff data:
numfuns = numel(f.funs);
n = cell(numfuns, 1);
c = get(f, 'coeffs');
% [TODO]: GET('coeffs') should not be used. We proably need chebpolyplotData()
% methods, similar to those used by the PLOT() funtion.
if ( ~iscell(c) )
    c = {c};
end
c = cellfun(@abs, c, 'UniformOutput', false);
for k = 1:numfuns
    n{k} = (size(c{k}, 1)-1:-1:0).';
end
numcols = size(c{1}, 2);

% Shape it:
data = reshape([n c]', 1, 2*numfuns);

% Plot the coeffs:
h1 = semilogy(data{:}, varargin{:});
hold on

% Get data for epslevel plot:
v = get(f, 'vscale-local');
e = get(f, 'epslevel-local');
ve = v.*repmat(e, 1, numcols);

% Reshape data for epslevel plot:
n = cellfun(@(x) x([end ; 1]), n, 'UniformOutput', false);
n = reshape(repmat(n', numcols, 1), numcols*numel(n), 1);
ve = reshape(ve, numfuns*numcols, 1);
ve = mat2cell(repmat(ve, 1, 2), ones(numfuns*numcols, 1), 2);
data = reshape([n ve]', 1, 2*numfuns*numcols);

% Plot the epslevels:
h2 = semilogy(data{:}, varargin{:});
hold off
for k = 1:numel(h2)
    c = get(h1(k), 'color');
    set(h2(k), 'linestyle', ':', 'linewidth', 1, 'marker', 'none', 'color', c);
end

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output if one was requested:
if ( nargout > 0 )
    varargout = {h1, h2};
end
