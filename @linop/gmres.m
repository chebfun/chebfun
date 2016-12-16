function varargout = gmres(A, varargin)
%GMRES   Iterative solution of a linear system.
%   U = GMRES(A,F) solves the system A*U = F for CHEBFUNs U and F and LINOP A.
%
%   More calling options are available; see CHEBFUN/GMRES for details.
%
% Example:
%   % To solve a simple Volterra integral equation:
%   d = domain(-1,1);
%   f = chebfun('exp(-4*x.^2)', d);
%   A = cumsum(d) + 20;
%   u = gmres(A, f, Inf, 1e-14);
%
% See also CHEBFUN/GMRES, GMRES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

op = @(u) A*u;
[varargout{1:nargout}] = gmres(op, varargin{:});

end   

