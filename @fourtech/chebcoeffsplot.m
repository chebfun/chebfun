function varargout = chebcoeffsplot(f, varargin)
%CHEBCOEFFSPLOT   Plot Chebyshev coefficients of a FOURTECH.
%
% See also COEFFSPLOT, FOURCOEFFSPLOT, PLOT

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Update documentation.

f = four2cheb(f);
[varargout{1:nargout}] = chebcoeffsplot(f, varargin{:});

end
