function varargout = cos(varargin)
%COS   Cosine of a SPHEREFUN.
%   COS(F) returns the cosine of F.
%
% See also SPHEREFUN/COSH and SPHEREFUN/SIN

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = cos@separableApprox(varargin{:});

end
