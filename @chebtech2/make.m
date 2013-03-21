function f = make(varargin)
%MAKE   Constructor shortcut for CHEBTECH2 objects.
%
%   For CHEBTECH2 methods implemented at the CHEBTECH level, it is not possible to
%   call the class constructor file corresponding to a CHEBTECH object directly.
%   F = MAKE(VARARGIN) allows us to get around this and construct a CHEBTECH2.

f = chebtech2(varargin{:}); 

end
