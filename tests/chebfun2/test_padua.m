function pass = test_padua( pref )

if ( nargin == 0 )
    pref = chebfunpref();
end

tol = 1e-13;
seedRNG(0)

%% Test ordering of PADUAPTS is consistent with PADPTS():

x = paduapts(2);
pdpts = [   1    0.5
            1   -1
            0    1
            0   -0.5
           -1    0.5
           -1   -1   ];
err(1) = norm(x - pdpts);

%% Test construction from PADPUAPTS():

% Small example on random data:
dom = [-2, 2, -2, 7];
F = chebfun2(rand(4), dom);
x = paduapts(6, dom);
f = feval(F, x(:,1), x(:,2));
G = chebfun2(f, dom, 'padua');
err(2) = norm(F - G);

% Larger example:
dom = [-1.3, 1, -1, 1.5];
FF = @(x, y) sin((x+.3).*(y-.2) + x.^2 + .4);
x = paduapts(150, dom);
f = feval(FF, x(:,1), x(:,2));
G = chebfun2(f, dom, 'padua');
err(3) = norm(FF(0,0) - G(0,0));

%%

pass = err < tol;

end