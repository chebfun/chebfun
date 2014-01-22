function out = roots(f, varargin)
%ROOTS   Roots of a DELTAFUN F.
%   ROOTS(F) returns the real roots of the FUNPART of the DELTAFUN F in the 
%   domain of F.
%
%   ROOTS(F, PROP1, VAL1, PROP2, VAL2, ...) modifies the default ROOTS
%   properties. The PROPs (strings) and VALs may be any of the following:
%
% See chebfun/roots for all the options.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with tne empty case:
if ( isempty(f) )
    out = [];
    return
end

out = roots(f.funPart, varargin{:} );
end