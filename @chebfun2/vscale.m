function varargout = vscale(varargin)
%VSCALE   Vertical scale of a CHEBFUN2.
%
% VSCL = VSCALE(F) returns the vertical scale of a CHEBFUN2 as determined
% by evaluating on a coarse tensor-product grid.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = vscale@separableApprox(varargin{:});

end
