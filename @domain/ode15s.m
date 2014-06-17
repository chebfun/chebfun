function varargout = ode15s(varargin)
%ODE15s   Solve stiff differential equations and DAEs. Output a CHEBFUN.
%   
% This syntax is deprecated. Please use chebfun.ode15s(...) instead.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:DOMAIN:ode15s:deprecated', ...
    ['Usage of ODE15S via the @DOMAIN class is deprecated and may be ', ...
    'removed from future releases. Please use chebfun.ode15s(...) instead.']);

varargin = domain.toDouble(varargin{:});

varargout{1:nargout} = chebfun.ode15s(varargin{:});

end
