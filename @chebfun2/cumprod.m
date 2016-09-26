function varargout = cumprod(varargin)
%CUMPROD  Indefinite product integral of a CHEBFUN2.
%   G = CUMPROD(F) returns the CHEBFUN2 G = exp( cumsum(log(F)) )
%
%   G = CUMPROD(F, DIM) returns the CHEBFUN2 G = exp( cumsum(log(F), DIM) )
%
% See also CUMSUM, SUM, PROD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = cumprod@separableApprox(varargin{:});

end
