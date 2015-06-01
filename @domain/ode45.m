function varargout = ode45(varargin)
%ODE45   Solve non-stiff ODEs. Output a CHEBFUN.
%   
% This syntax is deprecated. Please use chebfun.ode45(...) instead.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:DOMAIN:ode45:deprecated', ...
    ['Usage of ODE45 via the @DOMAIN class is deprecated and may be ', ...
    'removed from future releases. Please use chebfun.ode45(...) instead.']);

varargin = domain.toDouble(varargin{:});

varargout{1:nargout} = chebfun.ode45(varargin{:});

end
