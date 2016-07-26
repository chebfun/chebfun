function varargout = squeeze(varargin)
%SQUEEZE   Squeeze a DISKFUN to one variable, if possible.
%   G = squeeze(F) returns a DISKFUN if F depends on LAMBDA and THETA in
%   spherical coordinates. If F depends only on the lambda-variable a row
%   CHEBFUN is returned and if it depends on just the theta-variable a
%   column CHEBFUN is returned.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = squeeze@separableApprox(varargin{:});
end
