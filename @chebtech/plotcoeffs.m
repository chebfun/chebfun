function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display Chebyshev coefficients graphically.
%   PLOTCOEFFS(F) plots the Chebyshev coefficients of a CHEBTECH F on a semilogy
%   scale. A horizontal line at the EPSLEVEL of F is also plotted. If F is an
%   array-valued CHEBTECH then a curve is plotted for each component (column) of
%   F.
%
%   PLOTCOEFFS(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale. If S
%   contains a string 'NOEPSLEVEL', the EPSLEVEL is not plotted.
%
%   H = PLOTCOEFFS(F) returns a column vector of handles to lineseries objects.
%   The final entry is that of the EPSLEVEL plot.
%
%   Note: to make the PLOTCOEFFS easier to read, zero coefficients have a small
%   value added to them (typically F.EPSLEVEL) so that the curve displayed is
%   continuous.
%
% See also CHEBCOEFFS, PLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  The undocumented featrue plotcoeffs(f, 'barplot') shows a different kind of
%  coeffs plot, which can be more attractive in some situations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

% Set defaults:
loglogPlot = false;
plotEpsLevel = true;
doBar = false;

% Copy input arguments:
args = varargin;

% Check inputs for 'loglog' or 'noepslevel'.
j = 1;
while ( j <= length(args) )
    if ( strcmpi(args{j}, 'loglog') )
        loglogPlot = true; 
        args(j) = [];
    elseif ( strcmpi(args{j}, 'noepslevel') )
        plotEpsLevel = false; 
        args(j) = [];
    elseif ( strcmpi(args{j}, 'barplot') )
        doBar = true;
        args(j) = [];
    else
        j = j + 1;
    end
end

% Store the hold state of the current axis:
holdState = ishold;

% The coefficients:
absCoeffs = abs(f.coeffs);

% Add a tiny amount to zeros to make plots look nicer:
if ( f.vscale > 0 )
    if ( doBar )
        absCoeffs(absCoeffs < min(f.epslevel.*f.vscale)/100) = 0;
    else
        % Min of epslevel*vscale and the minimum non-zero coefficient:
        absCoeffs(~absCoeffs) = min( min(f.epslevel.*f.vscale), ...
                                 min(absCoeffs(logical(absCoeffs))) );                             
    end
else
    % Add epslevel for zero CHEBTECHs:
    absCoeffs = absCoeffs + f.epslevel;
end

% Get the size:
[n, m] = size(absCoeffs);

xx = 0:1:n-1;
yy = absCoeffs;
if ( any(doBar) )
    [xx, yy] = padData(xx,yy);
end

% Plot the coeffs:
h = semilogy(xx, yy, args{:}); 
hold on

if ( plotEpsLevel )
    % Plot the epslevel:
    h2 = semilogy([0 n-1], repmat(f.vscale.*f.epslevel, 2, 1), args{:});
    for k = 1:m
        c = get(h(k), 'color');
        set(h2(k), 'linestyle', ':', 'linewidth', 1, 'marker', 'none', 'color', c);
    end
    set(h2, 'handlevis', 'off');
else
    h2 = plot([]);
    set(h2, 'handlevis', 'off');
end

% Do a loglog plot:
if ( loglogPlot )
    set(gca, 'XScale', 'log')
end

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Adjust xLim:
xLim = get(gca, 'xlim');
set(gca, 'xLim', [min(xLim(1), 0), max(xLim(2), n)])

% Add title and labels
title(gca, 'Chebyshev coefficients')
xlabel(gca, 'Degree of Chebyshev polynomial')
ylabel(gca, 'Magnitude of coefficient')

% By default, set grid on
grid(gca, 'on')

% Give an output if one was requested:
if ( nargout > 0 )
    varargout{1} = h;
    varargout{2} = h2;
end

end

function [xx, yy] = padData(x, y)
%PADDATA   Pad the x and y data to make a bar plot:
xx = [x+.5 ; x-.5 ; x-.5];
xx(xx<0) = 0;
xx = xx(:);

[n, m] = size(y);
nans = NaN(n, 1);
yy = zeros(3*n, m);
for k = 1:size(y,2)
    yk = y(:,k);
    yk = [yk yk nans].';
    yy(:,k) = yk(:);
end
end
