function varargout = imag(varargin)
%IMAG   Imaginary part of a CHEBFUN2.
%   IMAG(F) returns the imaginary part of a CHEBFUN2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = imag@separableApprox(varargin{:});

end
