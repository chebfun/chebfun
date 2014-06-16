% Test file for JACPOLY.

function pass = test_jacpoly(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% TODO: We simply test against some values from V4. This should be improved.

tol = 1e-14;

J = jacpoly(0:3, .1, -.3);
Jx = feval(J, .1).';
v4 = [     1.000000000000000
           0.290000000000000
          -0.413700000000000
          -0.254186000000000];
err(1) = norm(Jx - v4, inf);

J = jacpoly(0:3, .1, -.3, [0, 3]);
Jx = feval(J, .1).';
v4 = [     1.000000000000000
          -0.640000000000000
           0.442244444444444
          -0.2715069629629630];

err(2) = norm(Jx - v4, inf);

pass = err < tol;

end