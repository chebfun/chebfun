function varargout = power(varargin)
%.^	     Pointwise power of a DISKFUN. 
%
% F.^G returns a DISKFUN F to the scalar power G, a scalar F to the
% DISKFUN power G, or a DISKFUN F to the DISKFUN power G.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = power@separableApprox(varargin{:});

end