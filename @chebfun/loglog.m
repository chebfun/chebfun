function varargout = loglog(varargin)
%LOGLOG  log-log scale plot of a CHEBFUN.
%   LOGLOG(...) is the same as PLOT(...), except a logarithmic (base 10) scale
%   is used for the X- and Y-axes.
%
% See also PLOT, SEMILOGX, SEMILOGY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Standard chebfun/plot:
h = plot(varargin{:});

% Set the XScale and YScale to be logarithmic:
set(gca, 'XScale', 'log', 'YScale', 'log');

% Output handle if requested:
if ( nargout > 0 )
    varargout = {h};
end

end
