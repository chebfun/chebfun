function varargout = squeeze(varargin)
%SQUEEZE   Squeeze a CHEBFUN2 to one variable, if possible.
%   G = squeeze(F) returns a CHEBFUN2 if F depends on x and y. If F depends only
%   on the x-variable a row CHEBFUN is returned and if it depends on just the
%   y-variable a column CHEBFUN is returned.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = squeeze@separableApprox(varargin{:});

end
