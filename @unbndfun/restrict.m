function g = restrict(f, s, pref)
%RESTRICT Restrict an UNBNDFUN to a subinterval.
%   RESCTRICT(F, S) returns a UNBNDFUN that is restricted to the subinterval
%   [S(1), S(2)] of F.domain.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESCTRICT(F, S) returns
%   an array of FUN objects, where the cells contain F restricted to each of 
%   the subintervals defined by S. If there is only one FUN to be returned,
%   that is, length(S) == 2, then the FUN object g is returned. This 
%   facilitates the use of the result by other functions, e.g. plot etc.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Check if subint is actually a subinterval:
if ( s(1) < f.domain(1) || s(end) > f.domain(2) || any(diff(s) <= 0) )
    error('UNBNDFUN:restrict:badinterval', 'Not a valid interval.')
elseif ( numel(s) == 2 && all(s == f.domain) )
    % Nothing to do here!
    return
end

% grab some preferneces:
if ( nargin < 3 ) 
    pref = fun.pref;
end
    
% Grab the scales:
hscale = get(f, 'hscale');
vscale = get(f, 'vscale');

% Loop over each of the new subintervals and make a new FUN:
g = cell(1, numel(s) - 1);
for k = 1:(numel(s) - 1)
    g{k} = fun.constructor(@(x) feval(f, x), s(k:k+1), hscale, vscale, pref);
end

end