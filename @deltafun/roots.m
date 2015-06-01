function out = roots(f, varargin)
%ROOTS   Roots of a DELTAFUN F.
%   ROOTS(F) returns the real roots of the FUNPART of the DELTAFUN F in the
%   domain of F.
%
%   ROOTS(F, PROP1, VAL1, PROP2, VAL2, ...) modifies the default ROOTS
%   properties. See the documentation in FUNPART/ROOTS() for more details.
%
% See chebfun/roots for all the options.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

out = roots(f.funPart, varargin{:});

end
