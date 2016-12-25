function varargout = rdivide(varargin)
%./   Pointwise right divide of DISKFUN objects.
%   F./G returns (1/G)*F if F is a DISKFUN and G is a double.
%
%   F./G returns F/G if F is a double and G is a DISKFUN. This does not 
%   work if G becomes numerically close to zero.
% 
%   F./G is the same as the command rdivide(F,G)
%
% See also DISKFUN/LDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = rdivide@separableApprox(varargin{:});

end