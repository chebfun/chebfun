function varargout = ode113(varargin)
%ODE113   Solve stiff differential equations and DAEs. Output a CHEBFUN.
%   
% This syntax is depricated. Please use chebfun.ode113(...) instead.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:ode113:depricated', ...
    ['Usage of ODE113 via the @DOMAIN class is depreicated and may be ', ...
    'removed from future releases. Please use chebfun.ode113(...) instead.']);

varargin = domain.toDouble(varargin{:});

varargout{1:nargout} = chebfun.ode113(varargin{:});

end
