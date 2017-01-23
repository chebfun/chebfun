function varargout = feval(varargin)
%FEVAL   Evaluate a CHEBFUN2 at one or more points.
%   FEVAL(F, X, Y) evaluates the CHEBFUN2 F at the point(s) in (X, Y), 
%   where X and Y are doubles.
%
%   FEVAL(F, X) evaluates the CHEBFUN2 F at the point(s) in X, where X is a
%   double, interpreted as F(real(X), imag(X)).
%
% See also SUBSREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = feval@separableApprox(varargin{:});

end