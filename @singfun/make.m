function f = make(varargin)
%MAKE   Constructor shortcut for SINGFUN objects.
%   For SINGFUN methods implemented at the ONEFUN level, it is not possible to 
%   call the class constructor file corresponding to a ONEFUN object directly. 
%   F = MAKE(VARARGIN) allows us to get around this and construct a SINGFUN.
%
% See also SINGFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Developer note: This is essentially a "factory method" in the sense of
% _Design Patterns_ by Gamma, Helm, Johnson, and Vlissides.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = singfun(varargin{:}); 

end
