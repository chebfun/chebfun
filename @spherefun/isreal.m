function varargout = isreal(varargin)
%ISREAL   Real-valued SPHEREFUN2 test.
%   ISREAL(F) returns logical true if F does not have an imaginary part and
%   false otherwise.
%  
%   ~ISREAL(F) detects SPHEREFUN2 object that have an imaginary part even if it i
%   all zero.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = isreal@separableApprox(varargin{:});
end
