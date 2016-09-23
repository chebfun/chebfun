function varargout = isreal(varargin)
%ISREAL   Real-valued CHEBFUN2 test.
%   ISREAL(F) returns logical true if F does not have an imaginary part and
%   false otherwise.
%
%   ~ISREAL(F) detects CHEBFUN2 object that have an imaginary part even if it is
%   all zero.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = isreal@separableApprox(varargin{:});

end
