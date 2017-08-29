function varargout = vscale(varargin)
%VSCALE   Vertical scale of a SPHEREFUN.
% 
% VSCL = VSCALE(F) returns the vertial scale of a SPHEREFUN as determined
% by evaluating on a coarse tensor-product grid. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = vscale@separableApprox(varargin{:});

end
