function varargout = plotcoeffs2(varargin)
%PLOTCOEFFS2   Display bivariate coefficients graphically.
%   PLOTCOEFFS2(F) plots the bivariate coefficients in a stem3 plot with a
%   semilogy scale.
%
%   H = PLOTCOEFFS2(F) returns a handle H to the figure.
%
% See also PLOTCOEFFS, CHEBCOEFFS2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = plotcoeffs2@separableApprox(varargin{:});

end
