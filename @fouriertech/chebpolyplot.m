function varargout = chebpolyplot(f, varargin)
%CHEBPOLYPLOT   Display Fourier coefficients graphically.
%
%   CHEBPOLYPLOT(F) plots the Fourier coefficients of a FOURIERTECH F on a
%   semilogy scale. A horizontal line at the EPSLEVEL of F is also plotted. If F
%   is an array-valued FOURIERTECH then a curve is plotted for each component
%   (column) of F.
%
%   CHEBPOLYPLOT(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale. If S
%   contains a string 'NOEPSLEVEL', the EPSLEVEL is not plotted.
%
%   H = CHEBPOLYPLOT(F) returns a column vector of handles to lineseries
%   objects. The final entry is that of the EPSLEVEL plot.
%
%   Note: to make the CHEBPOLYPLOT easier to read, zero coefficients have a
%   small value added to them (typically F.EPSLEVEL) so that the curve displayed
%   is continuous.
%
% See also CHEBPOLY, PLOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

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
    % (Min of epslevel*vscale and the miniumum non-zero coefficient)
    absCoeffs(~absCoeffs) = min( f.epslevel*min(f.vscale), ...
                                 min(absCoeffs(logical(absCoeffs))) );
else
    % (add epslevel for zero FOURIERTECHs)
    absCoeffs = absCoeffs + f.epslevel;
end

% Get the size:
[n, m] = size(absCoeffs);

if ( plotEpsLevel )
    % Plot the coeffs AND the epslevel:
    h = semilogy((n-1)/2:-1:-(n-1)/2, absCoeffs, args{:});
    hold on
    h2 = semilogy([(n-1)/2 -(n-1)/2], repmat(f.vscale, 2, 1)*f.epslevel, args{:});
    h = [h ; h2];
    set(gca,'XDir', 'Reverse')
    for k = 1:m
        c = get(h(k), 'color');
        set(h(m+k), 'linestyle', ':', 'linewidth', 1, 'marker', 'none', 'color', c);
    end
else
    % Plot just the coefficients:
    h = semilogy((n-1)/2:-1:-(n-1)/2, absCoeffs, args{:});
end

% For constant functions, plot a dot:
if ( n == 1 )
    set(h, 'marker', 'o');
end

% Do a loglog plot:
if ( loglogPlot )
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
