function varargout = isempty(varargin)
%ISEMPTY   True for empty CHEBFUN2.
%   ISEMPTY(F) returns 1 if F is an empty CHEBFUN2 object and 0 otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = isempty@separableApprox(varargin{:});

end
