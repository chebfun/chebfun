function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display Trigonometric coefficients graphically.
%   PLOTCOEFFS(F) plots the Trigonometric coefficients of a TRIGTECH F on a
%   semilogy scale.  If F is an array-valued TRIGTECH then a curve is plotted
%   for each component (column) of F.
%
%   PLOTCOEFFS(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale.
%
%   H = PLOTCOEFFS(F) returns a column vector of handles to lineseries
%   objects.
%
%   Note: to make the COEFPLOT easier to read, zero coefficients have a small
%   value added to them (typically EPS*VSCALE(F)) so that the curve displayed
%   is continuous.
%
% See also PLOT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
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
domain = [-1, 1];

% Copy input arguments:
args = varargin;

% Check inputs for 'loglog'.
j = 1;
while ( j <= length(args) )
    if ( strcmpi(args{j}, 'loglog') )
        loglogPlot = true; 
        args(j) = [];
    elseif ( strcmpi(args{j}, 'domain') )
        domain = args{j+1};
        args(j:j+1) = [];
    else
        j = j + 1;
    end
end

% Store the hold state of the current axis:
holdState = ishold;

% The coefficients:
absCoeffs = abs(f.coeffs);

% Get the size:
[n, m] = size(absCoeffs);

% Need to handle odd/even cases separately
isEven = ~mod(n, 2);
if ( isEven )
    coeffIndex = -n/2:n/2-1;
else
    coeffIndex = -(n-1)/2:(n-1)/2;
end

% Add a tiny amount to zeros to make plots look nicer:
if ( vscale(f) > 0 )
    % (Min of eps*vscale and the minimum non-zero coefficient)
    absCoeffs(~absCoeffs) = min( min(eps*vscale(f)), ...
                                 min(absCoeffs(logical(absCoeffs))) );                             
else
    % (add eps for zero CHEBTECHs)
    absCoeffs = absCoeffs + eps;
end

if ( ~loglogPlot )
    % Plot the coefficients:
    normalizedWaveNumber = coeffIndex*(2*pi)/diff(domain);
    h = semilogy(normalizedWaveNumber, absCoeffs, args{:});
    if ( ~holdState )
        xlim([min(normalizedWaveNumber(1),-1) -min(normalizedWaveNumber(1),-1)]);
        xlim([min(normalizedWaveNumber(1),-1) -min(normalizedWaveNumber(1),-1)]);
    end
    % Set the string for the x-axis label.
    xlabelStr = 'Wave number';
else
    if ( isEven )
        % In this case the negative cofficients have an additional term
        % corresponding to the cos(N/2*x) coefficient. We will store
        % the constant mode coefficient in both vectors.
        cPos = absCoeffs(n/2+1:n,:);
        cNeg = absCoeffs(n/2+1:-1:1,:);
    else
        cPos = absCoeffs((n+1)/2:n,:);
        cNeg = absCoeffs((n+1)/2:-1:1,:);
    end
    coeffIndexPos = 1:size(cPos,1);
    coeffIndexNeg = 1:size(cNeg,1);
    waveNumber = [coeffIndexPos nan coeffIndexNeg];
    normalizedWaveNumber = waveNumber*(2*pi)/diff(domain);

    % Plot the coefficients for the positive and negative fourier modes
    % separately.
    h = loglog(normalizedWaveNumber, [cPos ; nan(1, m) ; cNeg], args{:});

    % Set the string for the x-axis label.  In this case we will be
    % plotting the absolute value of the wave number + 1 (since we can't
    % represent a wave number <= 0 on a logarithmic scale.
    xlabelStr = '|Normalized wave number|+1';
end

% For constant functions, plot a dot:
if ( n == 1 )
    set(h, 'marker', 'o');
end

% Add title and labels
title(gca, 'Fourier coefficients')
xlabel(gca, xlabelStr)
ylabel(gca, 'Magnitude of coefficient')

% By default, set grid on
grid(gca, 'on')

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output if one was requested:
if ( nargout > 0 )
    varargout{1} = h;
end

end
