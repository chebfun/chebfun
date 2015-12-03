function varargout = feval(varargin)
%FEVAL  Evaluate a CHEBFUN2 at one or more points.
%   FEVAL(F,X,Y) evaluates the CHEBFUN2 F and the point(s) in (X,Y), where X and
%   Y are doubles.
%
%   FEVAL(F,X) evaluates the CHEBFUN2 F along the complex valued CHEBFUN X and
%   returns  g(t) = F(real(X(t)),imag(X(t)))
%
%   FEVAL(F,X,Y) returns g(t) = F(X(t),Y(t)), where X and Y are real valued
%   CHEBFUN objects with the same domain.
%
% See also SUBSREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = feval@separableApprox(varargin{:});

end
