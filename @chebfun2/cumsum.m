function varargout = cumsum(varargin)
%CUMSUM   Indefinite integral of a CHEBFUN2.
%   F = CUMSUM(F) returns the indefinite integral of a CHEBFUN2 with respect to
%   one variable and hence, returns a chebfun. The integration is done by
%   default in the y-direction.
%
%   F = CUMSUM(F, DIM). If DIM = 2 integration is along the x-direction, if DIM
%   = 1 integration is along the y-direction.
%
% See also CUMSUM2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = cumsum@separableApprox(varargin{:});

end
