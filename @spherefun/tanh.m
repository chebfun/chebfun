function varargout = tanh(varargin)
%TANH   Hyperbolic tangent of a SPHEREFUN.
% 
%   TANH(F) returns the hyperbolic tangent of a SPHEREFUN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = tanh@separableApprox(varargin{:});
end
