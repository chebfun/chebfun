function pass = test_system3(pref)

if ( nargin == 0 )
    pref = cheboppref;
end

N = chebop(@(x,u,v) [diff(u,2) + sin(v) ; cos(u) + diff(v,2)]);
N.lbc = @(u,v) [u - 2 ; v - 1]; 
N.rbc = @(u,v) [u - 2 ; v + 1];
[u, v, info1] = solvebvp(N, [0 ; 0], pref);

pass = abs(feval(v, .2) - -0.371250985730553) < 1e-8;

end
