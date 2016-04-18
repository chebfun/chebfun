function varargout = cumprod(varargin)
%CUMPROD  Indefinite product integral of a SPHEREFUN. 
%   G = CUMPROD(F) returns the SPHEREFUN G = exp( cumsum(log(F)) )
% 
%   G = CUMPROD(F, DIM) returns the SPHEREFUN G = exp( cumsum(log(F), DIM) )
%
% See also CUMSUM, SUM, PROD.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = cumprod@separableApprox(varargin{:});
end
