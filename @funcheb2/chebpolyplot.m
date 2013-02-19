function varargout = chebpolyplot(f,varargin)
%CHEBPOLYPLOT    Display Chebyshev coefficients graphically.
%
%   CHEBPOLYPLOT(F) plots the Chebyshev coefficients of a FUNCHEB2 F on a
%   semilogy scale. A horizontal line at the epslevel of F is also plotted. If F
%   is a vectorised FUNCHEB2 then a curve is plotted for each component (column)
%   of F.
%
%   CHEBPOLYPLOT(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale. If S
%   contains a string 'NOEPSLEVEL' the epslevel is not plotted.
%
%   H = CHEBPOLYPLOT(F) returns a column vector of handles to lineseries
%   objects. The final entry is that of the epslevel plot.
%
% See also FUNCHEB2.CHEBPOLY, FUNCHEB2/PLOT

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Set defaults:
loglogplot = false;
plotepslevel = true;

% Copy input arguments:
args = varargin;

% Check inputs for 'loglog' or 'noepslevel'.
j = 1;
while ( j <= length(args) )
    if ( strcmpi(args{j},'loglog') )
        loglogplot = true; 
        args(j) = [];
    elseif ( strcmpi(args{j},'noepslevel') )
        plotepslevel = false; 
        args(j) = [];
    else
        j = j + 1;
    end
end

% Store the hold state of the current axis:
holdState = ishold;

% The coefficients:
absc = abs(f.coeffs);
% Add a tiny amount to zeros to make plots look nicer:
absc(~absc) = eps*max(absc(:));

% Get the size:
n = size(absc, 1);

if ( plotepslevel )
    % Plot the coeffs AND the epslevel:
    h = semilogy(n-1:-1:0, absc, [0 n-1], f.vscale*f.epslevel*[1, 1], args{:});
    set(h(end), 'linestyle', '--', 'linewidth', 1, 'color', 'r', ...
        'marker', 'none')
else
    % Plot just the coeffs:
    h = semilogy(n:-1:1, absc, args{:});
end

% Do a loglog plot:
if ( loglogplot )
    set(gca, 'XScale', 'log')
end

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output if one was requested:
if ( nargout > 0 )
    varargout{1} = h;
end

end

