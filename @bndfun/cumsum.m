function f = cumsum(f,m,varargin)
%CUMSUM   Indefinite integral of a BNDFUN.
%   CUMSUM(F) is the indefinite integral of the BNDFUN F on an interval [a,b],
%   with the constant of integration chosen so that F(a) = 0.
%
%   CUMSUM(F, M) will compute the Mth definite integral with the constant of
%   integration chosen so that each intermediate integral evaluates to 0 at x=a.
%   Thus CUMSUM(F, 2) is equivalent to CUMSUM(CUMSUM(F)).
%
%   CUMSUM(F, PREF) or CUMSUM(F, M,  PREF) uses options from the preference
%   structure PREF when building the output CHEBTECH.
%
% See also DIFF, SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

%%
% Trivial case of an empty BNDFUN:
if ( isempty(f) )
    return
end

% Parse inputs:
if ( nargin == 1 )
    m = 1;
    pref = bndfun.pref;
elseif ( nargin < 3 )
    if ( isstruct(m) )
        pref = m;
        m = 1;
    else
        pref = bndfun.pref();
    end
end

% Rescaling factor, (b-a)/2, to the mth power.
rescaleFactorm = (.5*diff(f.domain))^m;

% Compute the cumsum of all of f's onefuns, multiply by the rescaling factor,
% and assign to the onefun field of f.
f.onefun = cumsum(f.onefun,m,pref)*rescaleFactorm;

end
