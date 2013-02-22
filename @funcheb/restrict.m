function f = restrict(f, s, pref)
%RESTRICT Restrict a FUNCHEB to a subinterval.
%   RESCTRICT(F, S) returns a FUNCHEB that is restricted to the subinterval
%   [S(1),S(2)] of [-1, 1]. Note that that since FUNCHEB only live on [-1,1], a
%   linear change of variables is implicitly fapplied.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESCTRICT(F, S) returns
%   a multivalued FUNCHEB where each column of F is restricted to each of the
%   subintervals defined by S. If F is a vector-valued function, say [F1, F2],
%   then the restrict(F, S = [S1, S2, S3]) returns the vector-valued FUNCHEB
%   [restrict(F1,S). restrict(F2, S)].

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Check if subint is actually a subinterval
if ( s(1) < -1 || s(end) > 1 || any(diff(s) <= 0) )
    error('FUN:restrict:badinterval', 'Not a valid interval.')
elseif ( numel(s) == 2 && all(s == [-1, 1]) )
    % Nothing to do here!
    return
end

% Compute new values and coeffs:
n = length(f);
x = f.chebpts(n);                            % old grid
y = .5*[1-x, 1+x] * [s(1:end-1) ; s(2:end)]; % new grid
values = feval(f, y);                        % new values
coeffs = f.chebpoly(values);                 % new coeffs

% Append data to fun:
f.values = values;
f.coeffs = coeffs;

% Grab a tolerance:
if ( nargin < 3 )
    pref = f.pref(); 
end
pref.(class(f)).eps = max(pref.(class(f)).eps, f.epslevel);

% Simplify the output fun:
f = simplify(f, pref);

end

