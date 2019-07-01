function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display trigonometric coefficients graphically.
%   PLOTCOEFFS(F) plots the trigonometric coefficients of a TRIGTECH F on a
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
% See also TRIGCOEFFS, PLOT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
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
ms = NaN;

% Copy input arguments:
args = varargin;

% Check inputs for 'loglog' or 'domain' or 'markersize'
j = 1;
while ( j <= length(args) )
    if ( strcmpi(args{j}, 'loglog') )
        loglogPlot = true; 
        args(j) = [];
    elseif ( strcmpi(args{j}, 'domain') )
        domain = args{j+1};
        args(j:j+1) = [];
    elseif ( strcmpi(args{j}, 'markersize') )
        ms = args{j+1};
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

if ( vscale(f) == 0 )
    absCoeffs = absCoeffs + eps;
end

if isnan(ms)
    NN = length(coeffIndex);
    ms = 2.5 + 50/sqrt(NN+8);
end

linetype_specified = ( mod(length(args),2) == 1 );

if ( ~loglogPlot )
    % Plot the coefficients:
    normalizedWaveNumber = coeffIndex*(2*pi)/diff(domain);
    if linetype_specified
        h = semilogy(normalizedWaveNumber, absCoeffs, ...
                         args{1}, 'markersize', ms, args{2:end});
        warningFlag = strcmp(get(h,'Marker'),'none');
        if warningFlag
            diffVec = find(absCoeffs ~= 0);
            if ( length(diffVec) > 1 )
                if ( min(diff(diffVec)) >= 2 )
                warning('CHEBFUN:plotcoeffs', ['No lines will appear ', ...
                  'because of zero values. Use ''.'' or ''.-'' instead.']);
                end
            end
        end
    else
        h = semilogy(normalizedWaveNumber, absCoeffs, ...
                         '.', 'markersize', ms, args{:});
    end

    if ( ~holdState )
        xlim([min(normalizedWaveNumber(1),-1) -min(normalizedWaveNumber(1),-1)]);
        xlim([min(normalizedWaveNumber(1),-1) -min(normalizedWaveNumber(1),-1)]);
    end
    % Set the string for the x-axis label.
    xlabelStr = 'Wave number';
else
    if ( isEven )
        % In this case the negative coefficients have an additional term
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

    % Plot the coefficients for the positive and negative Fourier modes
    % separately.
%   h = loglog(normalizedWaveNumber, [cPos ; nan(1, m) ; cNeg], args{:});
    if linetype_specified
        h = semilogy(normalizedWaveNumber, [cPos ; nan(1, m) ; cNeg], ...
                         args{1}, 'markersize', ms, args{2:end});
    else
        h = semilogy(normalizedWaveNumber, [cPos ; nan(1, m) ; cNeg], ...
                         '.', 'markersize', ms, args{:});
    end

    % Set the string for the x-axis label.  In this case we will be
    % plotting the absolute value of the wave number + 1 (since we can't
    % represent a wave number <= 0 on a logarithmic scale.
    xlabelStr = '|Normalized wave number|+1';
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
