function f = make(varargin)
%MAKE   Constructor shortcut for CHEBTECH1 objects.
%   For CHEBTECH1 methods implemented at the CHEBTECH level, it is not possible
%   to call the class constructor file corresponding to a CHEBTECH object
%   directly. F = MAKE(VARARGIN) allows us to get around this and construct a
%   CHEBTECH1.

% See also CHEBTECH1.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = chebtech1(varargin{:}); 

end
