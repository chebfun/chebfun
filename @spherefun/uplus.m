function varargout = uplus(varargin)
%UPLUS   Unary plus for a SPHEREFUN2. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = uplus@separableApprox(varargin{:});
end
