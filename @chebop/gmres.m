function varargout = gmres(A, varargin)
%GMRES   Iterative solution of a linear system. 
%   U = GMRES(A,F) solves the system A*U = F for CHEBFUN U and F and linear 
%   CHEBOP A. If A is not linear, an error is returned.
%
%   More calling options are available; see chebfun/gmres for details.
%
% Example:
%   % To solve a simple Volterra integral equation:
%   d = [-1,1];
%   f = chebfun('exp(-4*x.^2)',d);
%   A = chebop(@(u) cumsum(u) + 20*u, d);
%   u = gmres(A,f,Inf,1e-14);
%
% See also CHEBFUN/GMRES, GMRES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~all(islinear(A)) )
    error('CHEBFUN:CHEBOP:gmres:nonlinear', ...
        'GMRES supports only linear CHEBOP instances.');
end

op = @(u) A*u;
[varargout{1:nargout}] = gmres(op, varargin{:});

end   

