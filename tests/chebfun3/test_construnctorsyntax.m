function pass = test_construnctorsyntax()
% This tests the Chebfun3 constructor for different syntax.

pass(1) = 1;
f = @(x,y,z) cos(x) + sin(x.*y.*z);  % simple function.
fstr = 'cos(x) + sin(x.*y.*z)';      % string version.

try
    % % Adaptive calls
    % Function handle
    chebfun3(f);
    
    % String
    chebfun3(fstr);
    
    % With domain.
    chebfun3(f, [-1 1 1 2 0 2]);
    
    % Split domain syntax
    chebfun3(f, [-1,0], [2,3]);
    
    % Function handle only in one variable.
    chebfun3(@(x,y,z) x);

    % Function handle only in two variables.
    chebfun3(@(x,y,z) x + cos(z));
catch
    pass(1) = 0 ;
end

% Test different ways to make Chebfun2 nonadaptive. 
dom = [-1 1 -1 1 2 5]; 
ff = @(x,y,z) sin(x.*y.*z);
[xx, yy, zz] = chebfun3.chebpts3(10, 15, 5, dom);
vals = ff(xx,yy,zz);
g = chebfun3(vals, dom);
[m, n, p] = length(g); % This won't work if 'simplify' is applied in the constructor.
pass(2) = m == 10 && n == 15 && p == 5;

% Construction with a given vector of lengths
g = chebfun3(ff, [10 15 5]);
[m, n, p] = length(g);
pass(3) = m == 10 && n == 15 && p == 5;

% Construction with a given vector of lengths on a given domain
g = chebfun3(ff, [10 15 5], [-2 1 -2 1 -2 1]);
[m, n, p] = length(g);
pass(4) = m == 10 && n == 15 && p == 5;
pass(5) = all(g.domain == [-2 1 -2 1 -2 1]);

% Construction with a given vector of lengths on a given domain with inputs
% in the reverse order:
g = chebfun3(ff, [-2 1 -2 1 -2 1], [10 15 5]);
[m, n, p] = length(g);
pass(6) = m == 10 && n == 15 && p == 5;
pass(7) = all(g.domain == [-2 1 -2 1 -2 1]);

% Construction with a specified rank
g = chebfun3(ff, 'rank', [2 3 4]);
[r1, r2, r3] = rank(g);
pass(8) = r1 == 2 && r2 == 3 && r3 == 4;
pass(9) = all(g.domain == [-1 1 -1 1 -1 1]);

% Construction with a specified rank on a given domain:
g = chebfun3(ff, 'rank', [2 3 4], [-2 1 -2 1 -2 1]);
[r1, r2, r3] = rank(g);
pass(10) = r1 == 2 && r2 == 3 && r3 == 4;
pass(11) = all(g.domain == [-2 1 -2 1 -2 1]);

% Construction with a specified rank on a given domain with inputs in the 
% reverse order:
g = chebfun3(ff, [-2 1 -2 1 -2 1], 'rank', [2 3 4]);
[r1, r2, r3] = rank(g);
pass(12) = r1 == 2 && r2 == 3 && r3 == 4;
pass(13) = all(g.domain == [-2 1 -2 1 -2 1]);

end