function varargout = imag(varargin)
%IMAG   Imaginary part of a SPHEREFUN.
%   IMAG(F) returns the imaginary part of a SPHEREFUN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = imag@separableApprox(varargin{:});
end
