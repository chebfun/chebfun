function f = make(varargin)
%MAKE   Constructor shortcut for UNBNDFUN objects.
%   For UNBNDFUN methods implemented at the FUN level, it is not possible to
%   call the class constructor file corresponding to a FUN object directly. 
%   F = MAKE(VARARGIN) allows us to get around this and construct a UNBNDFUN.

% See also UNBNDFUN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = unbndfun(varargin{:}); 

end
