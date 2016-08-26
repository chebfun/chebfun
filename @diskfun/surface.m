function varargout = surface(varargin)
%SURFACE  Plot surface of a DISKFUN.
%   See SURF for a complete description of the various forms of the input.
% 
% See also DISKFUN/SURF. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = surface@separableApprox(varargin{:});

end