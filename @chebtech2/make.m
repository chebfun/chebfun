function f = make(varargin)
%MAKE   Constructor shortcut for CHEBTECH2 objects.
%
%   For chebtech2 methods implented at the chebtech level, it is not possible to
%   call the class constructor file corresponding to a chebtech object directly.
%   F = MAKE(VARARGIN) allows us to get around this and construct a chebtech2.

f = chebtech2(varargin{:}); 

end
