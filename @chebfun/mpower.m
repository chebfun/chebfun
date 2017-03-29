function varargout = mpower(varargin)
%^   CHEBFUN power.
%   F^G is equivalent to F.^G or POWER(F, G);
%
% See also POWER.

% This is just a wrapper:
[varargout{1:nargout}] = power(varargin{:});

end