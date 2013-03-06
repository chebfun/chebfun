function g = restrict(f, s)
%RESTRICT Restrict a BNDFUN to a subinterval.
%   RESCTRICT(F, S) returns a BNDFUN that is restricted to the subinterval
%   [S(1), S(2)] of F.domain.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESCTRICT(F, S) returns
%   an array of FUNCHEB objects, where the cells contain F restricted to each of
%   the subintervals defined by S.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Check if subint is actually a subinterval:
if ( s(1) < f.domain(1) || s(end) > f.domain(2) || any(diff(s) <= 0) )
    error('BNDFUN:restrict:badinterval', 'Not a valid interval.')
elseif ( numel(s) == 2 && all(s == [-1, 1]) )
    % Nothing to do here!
    return
end

% Compute the breaks in [-1,1] space and restrict the onefun:
t = f.mapping.inv(s);
onefuns = restrict(f.onefun, t); % Returns an array of onefuns.

% Loop over each of the new subintervals and make a bndfun with new mapping:
g = cell(1, numel(s) - 1);
for k = 1:(numel(s) - 1)
    g{k} = bndfun(onefuns(k), s(k:k+1));
end

end

