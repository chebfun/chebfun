function varargout = times(varargin)
% .*    Pointwise multiplication for SPHEREFUN2 objects.
%
%   F.*G multiplies SPHEREFUN2 objects F and G. Alternatively F or G could be a
%   double.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = times@separableApprox(varargin{:});
end
