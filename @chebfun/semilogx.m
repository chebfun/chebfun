function varargout = semilogx(varargin)
%SEMILOGX   Semi-log scale plot of a CHEBFUN.
%   SEMILOGX(...) is the same as PLOT(...), except a logarithmic (base 10) scale
%   is used for the X-axis.
%
% See also PLOT, SEMILOGY, LOGLOG.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Standard chebfun/plot:
h = plot(varargin{:});

% Set the XScale to be logarithmic:
set(gca, 'XScale', 'log');

% Output handle if requested:
if ( nargout > 0 )
    varargout = {h};
end

end
