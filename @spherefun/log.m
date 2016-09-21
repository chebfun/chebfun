function varargout = log(varargin)
%LOG   Natural logarithm of a SPHEREFUN.
%   LOG(F) is the natural logarithm of F. This function returns an error 
%   if the function passes through or becomes numerically close to zero.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = log@separableApprox(varargin{:});

end
