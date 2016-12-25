function varargout = imag(varargin)
%IMAG   Imaginary part of a SPHEREFUN.
%   IMAG(F) returns the imaginary part of a SPHEREFUN.
%
%   Since only real-valued SPHEREFUNS are presently supported, this
%   function returns a zero spherefun.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = imag@separableApprox(varargin{:});

end
