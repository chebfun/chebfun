function pass = test_jumps_manual(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

% BVP with nonhomogeneous jump condition
% Toby Driscoll, June 2014
% Based on issue #665

tol = 1e1*pref.bvpTol;

%%
N = chebop(@(s, V) diff(V, 2), [-1, 1]);
N.lbc = @(V) V + 1;
N.rbc = @(V) V - 1;
N.bc = @(s,V) [jump(V, 0) - 1 ; jump(diff(V),0)];
y = N\0;

%%
err(1) = norm(N(y), inf);
jump1  = feval(y, 0, 'right') - feval(y, 0, 'left');
err(2) = abs(jump1 - 1);
jump2  = feval(diff(y), 0, 'right') - feval(diff(y), 0, 'left');
err(3) = abs(jump2);
err(4) = abs(y(-1) + 1);
err(5) = abs(y(1) - 1);

%%
pass = err < tol;

end
