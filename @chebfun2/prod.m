function varargout = prod(varargin)
%PROD     Product integral of a CHEBFUN2.
%   PROD(F) returns exp( sum(log(F)) )
%
%   PROD(F, DIM) returns the chebfun exp( sum(log(F), DIM) )
%
% See also CUMPROD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = prod@separableApprox(varargin{:});

end
