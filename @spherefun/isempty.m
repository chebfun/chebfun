function varargout = isempty(varargin)
%ISEMPTY   True for empty SPHEREFUN.
%   ISEMPTY(F) returns 1 if F is an empty SPHEREFUN object and 0 otherwise. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = isempty@separableApprox(varargin{:});

end
