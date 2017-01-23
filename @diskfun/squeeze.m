function varargout = squeeze(varargin)
%SQUEEZE   Squeeze a DISKFUN to one variable, if possible.
%   G = squeeze(F) returns a DISKFUN if F depends on THETA and R in
%   polar coordinates. If F depends only on the theta-variable, a row
%   CHEBFUN is returned. If it depends on just the r-variable, a
%   column CHEBFUN is returned.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = squeeze@separableApprox(varargin{:});

end