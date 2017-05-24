function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display Chebyshev coefficients graphically.
%   PLOTCOEFFS(F) plots the Chebyshev coefficients of a CHEBTECH F on a
%   semilogy scale.  If F is an array-valued CHEBTECH then a curve is plotted
%   for each component (column) of F.
%
%   PLOTCOEFFS(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale.
%
%   H = PLOTCOEFFS(F) returns a column vector of handles to lineseries objects.
%
% See also CHEBCOEFFS, PLOT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
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
doBar = false;
domain = [-1,1];
ms = NaN;

% Copy input arguments:
args = varargin;

% Check inputs for 'loglog' or 'barplot' or 'domain' or 'markersize'
j = 1;
while ( j <= length(args) )
    if ( strcmpi(args{j}, 'loglog') )
        loglogPlot = true; 
        args(j) = [];
    elseif ( strcmpi(args{j}, 'barplot') )
        doBar = true;
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

% The coefficients and vertical scale:
absCoeffs = abs(f.coeffs);
vscl = vscale(f);

% Add a tiny amount to zeros to make plots look nicer:
if ( vscl > 0 )
    if ( doBar )
        absCoeffs(absCoeffs < min(eps*vscl)/100) = 0;
    end
else
    % Add eps for zero CHEBTECHs:
    absCoeffs = absCoeffs + eps;
end

% Get the size:
[n, m] = size(absCoeffs);

xx = 0:1:n-1;
yy = absCoeffs;
if ( any(doBar) )
    [xx, yy] = padData(xx,yy);
end

% Plot the coefficients:
if isnan(ms)
    NN = length(xx);
    ms = 2.5 + 50/sqrt(NN+8);
end
linetype_specified = ( mod(length(args),2) == 1 );
if linetype_specified
    h = semilogy(xx, yy, args{1}, 'markersize', ms, args{2:end});
    warningFlag = strcmp(get(h,'Marker'),'none');
    if warningFlag
        diffVec = find(yy ~= 0);
        if ( length(diffVec) > 1 )
            if ( min(diff(diffVec)) >= 2 )
                warning('CHEBFUN:plotcoeffs', ['No lines will appear ', ...
                  'because of zero values. Use ''.'' or ''.-'' instead.']);
            end
        end
    end
else
    h = semilogy(xx, yy, '.', 'markersize', ms, args{:}); 
end
hold on

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
