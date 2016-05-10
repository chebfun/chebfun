function varargout = pivotplot( f, varargin )
%PIVOTPLOT   Semilogy plot of pivot values.
%   PIVOTPLOT( F ) semilogy plot related to the pivots values taken during
%   the construction of the SEPARABLEAPPROX F.
%
%   H = PIVOTPLOT( F ) returns a handle H to the figure.
%
%   PIVOTPLOT( F, S ) allows further plotting options, such as linestyle,
%   linecolor, etc. If S contains a string 'LOGLOG', the pseudo sig will be
%   displayed on a log-log scale.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) ) 
    error('CHEBFUN:SEPARABLEAPPROX:pivotplot:empty', ...
        'Empty SEPARABLEAPPROX has no pivots to plot');
end 

% Parse input arguments:
scaleTypeLog = false;
if ( nargin > 1 )
    for j = 1:length(varargin)
        if ( strcmpi(varargin{j}, 'loglog') )
            scaleTypeLog = true;
            varargin(j) = [];
            break
        end
    end
end

% Plot a semilogy or loglog: 
plotData = abs(pivots(f));
holdState = ishold;
if ( ~scaleTypeLog )
    h = semilogy( plotData, varargin{:} );            % semilogy plot
else
    h = loglog( plotData, varargin{:} );              % loglog plot
end

if ( ~holdState )
    hold off
end

% Output handle if appropriate:
if ( nargout > 0 )
    varargout = { h };
end

end
