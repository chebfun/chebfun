function f = make(varargin)
%MAKE   Constructor shortcut for FOURIERTECH objects.
%   For FOURIERTECH methods implemented at the FOURIERTECH level, it is not possible
%   to call the class constructor file corresponding to a FOURIERTECH object
%   directly. F = MAKE(VARARGIN) allows us to get around this and construct a
%   FOURIERTECH.
%
% See also FOURIERTECH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Developer note: This is essentially a "factory method" in the sense of
% _Design Patterns_ by Gamma, Helm, Johnson, and Vlissides.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fouriertech(varargin{:}); 

end
