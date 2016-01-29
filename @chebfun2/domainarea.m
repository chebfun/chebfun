function varargout = domainarea(varargin)
%DOMAINAREA    Area of the domain of f
%
%   DOMAINAREA(F) returns the area of the topological domain of f.
%

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = domainarea@separableApprox(varargin{:});

end
