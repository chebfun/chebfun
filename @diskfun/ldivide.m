function varargout = ldivide(varargin)
%.\   Pointwise DISKFUN left array divide.
%   F.\G if G is a DISKFUN and F is a double this returns (1/F)*G
%
%   F.\G if G is a double and F is a DISKFUN this returns G\F, but this
%   does not work if F becomes numerically close to zero.
%
%   F.\G we do not allow F and G to both be DISKFUN objects.
% 
%   F.\G is the same as the command ldivide(F,G)

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = ldivide@separableApprox(varargin{:});

end
