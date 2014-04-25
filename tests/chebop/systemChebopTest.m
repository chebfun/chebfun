function pass = systemChebopTest

% Setup
dom = [-1 1];
p = cheboppref;
p.plotting = 'off';
p.display = 'off';
p.damped = 1;

% cheboppref('display','iter')
% cheboppref('plotting','on')
N = chebop(@(x,u,v) [diff(u,2) + sin(v) ; cos(u) + diff(v,2)]);
N.lbc = @(u,v) [u - 2 ; v - 1]; 
N.rbc = @(u,v) [u - 2 ; v + 1];
[uv, info1] = solvebvp(N, [0 ; 0], p);

pass = abs(feval(uv{2}, .2) - -0.371250985730553) < 1e-8;

end