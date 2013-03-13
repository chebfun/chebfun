function f = make(varargin)
%MAKE   Constructor shortcut for FUNCHEB2 objects.
%
%   For funcheb2 methods implented at the funcheb level, it is not possible to
%   call the class constructor file corresponding to a funcheb object directly.
%   F = MAKE(VARARGIN) allows us to get around this, and construct a funcheb2.

f = funcheb2(varargin{:}); 

end