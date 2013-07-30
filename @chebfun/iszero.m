function out = iszero(f, tol)

% Grab a preference:
if ( nargin < 2 )
    vs = get(f, 'vscale'); vs = max([vs{:}]);
    el = max(get(f, 'epslevel'));
    tol = el*vs;
end

% Impulses:
out = all(f.impulses(:,1) < tol);

% Loop over each of the funs:
k = 0;
while ( k < numel(f.funs) && out )
    k = k + 1;
    out = iszero(f.funs{k});
end

end