function pass = test_feval2(pref)

pass = [];
if ( nargin == 0 )
    pref = chebfunpref();
end

%% Check chebop application
op  = @(x,u) diff(u,2) + x.*u;
bc  = @(x,u) [ u(-1); u(1) ];
A = chebop(op, [-1,1], bc);
[ev, ew] = eigs(A);
pass(1) = doesNotCrash(@() A*chebfun(ev));
pass(2) = doesNotCrash(@() A*quasimatrix(ev));
pass(3) = doesNotCrash(@() A*ev);

end


%% A function to test whether an expression crashes or not.
% Taken from tests/chebfun/test_plot.m
function pass = doesNotCrash(fn)

try
    fn();
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end

end
