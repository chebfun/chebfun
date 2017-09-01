function varargout = imag(varargin)
%IMAG   Imaginary part of a DISKFUN.
%   IMAG(F) returns the imaginary part of a DISKFUN.
%
%   Since only real-valued DISKFUNS are presently supported, this
%   function returns a zero diskfun.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = imag@separableApprox(varargin{:});

end
