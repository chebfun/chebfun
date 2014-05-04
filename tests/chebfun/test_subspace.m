% Test file for @chebfun/subspace.m.

function pass = test_subspace(pref)

if ( nargin == 0 )
   pref = chebfunpref();
end

% Orthonormal array-valued chebfun:
A = chebfun(@(t) [1/sqrt(2)+0*t cos(t) sin(2*t) sin(3*t)]/sqrt(pi), ...
    [0 1 2*pi], pref);

f = chebfun(@(t) sin(10*t)/sqrt(pi), [0 2 2*pi], pref);
alpha = [1e-10 pi/5 pi/2-1e-10];

for k = 1:length(alpha)
    B = cos(alpha(k))*A(:,k) + sin(alpha(k))*f;
    angle = subspace(A, B);
    pass(k) = abs(angle - alpha(k)) < 1e3*epslevel(B);
    angle = subspace(B, A);
    pass(k) = pass(k) && (abs(angle - alpha(k)) < 1e3*epslevel(B));
end

% Check subspaces with multiple columns.
pass(4) = abs(subspace(A(:,1:2), A(:,3:4)) - pi/2) < 10*vscale(A)*epslevel(A);
pass(5) = abs(subspace(A(:,1:3), A(:,4)) - pi/2) < 10*vscale(A)*epslevel(A);

% Check operation for row chebfuns.
At = A.';
pass(6) = abs(subspace(At(1:2,:), At(3:4,:)) - pi/2) < ...
    10*vscale(At)*epslevel(At);

end
