function varargout = cos(varargin)
%COS   Cosine of a DISKFUN.
%   COS(F) returns the cosine of F.
%
% See also DISKFUN/COSH and DISKFUN/SIN

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = cos@separableApprox(varargin{:});

end
