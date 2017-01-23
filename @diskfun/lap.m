function varargout = lap(varargin)
%LAP   Laplacian of a DISKFUN.
%   L = LAP(F) returns a DISKFUN representing the Laplacian of F. 
%
%   This is shorthand for LAPLACIAN(F).
%
% See also DISKFUN/LAPLACIAN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = laplacian(varargin{:});

end
