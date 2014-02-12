function varargout = ode45(varargin)
%ODE45   Solve non-stiff ODEs. Output a CHEBFUN.
%   
% This syntax is depricated. Please use chebfun.ode45(...) instead.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:ode113:depricated', ...
    ['Usage of ODE45 via the @DOMAIN class is depreicated and may be ', ...
    'removed from future releases. Please use chebfun.ode45(...) instead.']);

varargin = domain.toDouble(varargin{:});

varargout{1:nargout} = chebfun.ode45(varargin{:});

end
