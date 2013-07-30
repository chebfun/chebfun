function varargout = semilogy(varargin)
%SEMILOGY Semi-log scale plot of a chebfun.
%   SEMILOGY(...) is the same as PLOT(...), except a logarithmic (base 10) scale
%   is used for the Y-axis.
%
% See also PLOT, SEMILOGX, LOGLOG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Standard chebfun/plot:
h = plot(varargin{:});

% Set the YScale to be logarithmic:
set(gca, 'YScale', 'log');

% Output handle if requested:
if ( nargout > 0 )
    varargout = {h};
end

end