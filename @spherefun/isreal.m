function varargout = isreal(varargin)
%ISREAL   Real-valued SPHEREFUN test.
%   ISREAL(F) returns logical true if F does not have an imaginary part and
%   false otherwise.
%  
%   ~ISREAL(F) detects SPHEREFUN object that have an imaginary part even if it i
%   all zero.
%
%   Since only real-valued SPHEREFUNS are presently supported, this
%   function always returns a logical true.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = isreal@separableApprox(varargin{:});

end
