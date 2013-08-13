function out = iszero(f, tol)

% Grab a tolerance:
if ( nargin < 2 )
    vs = vscale(f);
    el = epslevel(f);
    tol = el*vs;
end

% Impulses:
out = all(f.impulses(:,1) < tol);

% Loop over each of the FUNs:
k = 0;
while ( k < numel(f.funs) && out )
    k = k + 1;
    out = iszero(f.funs{k});
end

end