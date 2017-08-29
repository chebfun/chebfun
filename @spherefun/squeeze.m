function varargout = squeeze(varargin)
%SQUEEZE   Squeeze a SPHEREFUN to one variable, if possible.
%   G = squeeze(F) returns a SPHEREFUN if F depends on LAMBDA and THETA in
%   spherical coordinates. If F depends only on the lambda-variable a row
%   CHEBFUN is returned and if it depends on just the theta-variable a
%   column CHEBFUN is returned.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = squeeze@separableApprox(varargin{:});
end
