function varargout = log(varargin)
%LOG   Natural logarithm of a CHEBFUN2.
%   LOG(F) is the natural logarithm of F. This function does not work if the
%   function passes through or becomes numerically close to zero.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = log@separableApprox(varargin{:});

end
