function f = make(varargin)
%MAKE   Constructor shortcut for TRIGTECH objects.
%   For TRIGTECH methods implemented at the TRIGTECH level, it is not possible
%   to call the class constructor file corresponding to a TRIGTECH object
%   directly. F = MAKE(VARARGIN) allows us to get around this and construct a
%   TRIGTECH.
%
% See also TRIGTECH.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Developer note: This is essentially a "factory method" in the sense of
% _Design Patterns_ by Gamma, Helm, Johnson, and Vlissides.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = trigtech(varargin{:}); 

end
