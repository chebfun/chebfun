function f = make(varargin)
%MAKE   Constructor shortcut for CHEBTECH2 objects.
%   For CHEBTECH2 methods implemented at the CHEBTECH level, it is not possible
%   to call the class constructor file corresponding to a CHEBTECH object
%   directly. F = MAKE(VARARGIN) allows us to get around this and construct a
%   CHEBTECH2.
%
% See also CHEBTECH2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Developer note: This is essentially a "factory method" in the sense of
% _Design Patterns_ by Gamma, Helm, Johnson, and Vlissides.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = chebtech2(varargin{:}); 

end
