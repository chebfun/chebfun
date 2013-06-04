function f = simplify(f, varargin)
%SIMPLIFY   Simplify the ONEFUN of the FUN object F. 
%   G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the FUN
%   object F, as defined by the SIMPLIFY method of the ONEFUN of F.
%
%   G = SIMPLIFY(F, PREF) does the same, but using some of the preferences in
%   the preference structure PREF.
%
% See also PREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Simplify the onefun of f.
f.onefun = simplify(f.onefun,varargin{:});
 
end
