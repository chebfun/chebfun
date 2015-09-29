% Test file for trigtech/diffmat.m

function pass = test_diffmat(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();
valsDiscretizationclass = trigcolloc();

%% 1st-order derivative:

%% Force n to be odd:
op = @(x) exp(cos(pi*x));
f = testclass.make(op, [], pref);
n = length(f);
if ( ~rem(n,2) )
    n = n+1;
    x = trigpts(n);
    v = op(x);
else
    x = trigpts(n);
    v = get(f, 'values');
end
D = valsDiscretizationclass.diffmat(n, 1);
df = D*v;
df_exact = @(x) -pi*sin(pi*x).*exp(cos(pi*x));
err = df_exact(x) - df;
pass(1) = ( norm(err, inf) < 1e3*vscale(f)*eps );
    


a = 10; b = 20;
f = testclass.make(@(x) cos(a*pi*sin(b*pi*x)), [], pref);
n = length(f);
if ( ~rem(n,2) )
    n = n+1;
    x = trigpts(n);
    v = op(x);
else
    x = trigpts(n);
    v = get(f, 'values');
end
D = valsDiscretizationclass.diffmat(n, 1);
df = D*v;
df_exact = @(x) -pi^2*a*b*cos(b*pi*x).*sin(a*pi*sin(b*pi*x));
err = df_exact(x) - df;
pass(2) = ( norm(err, inf) < 1e7*vscale(f)*eps );
    

%% Force n to be even:
op = @(x) exp(-50*x.^2);
f = testclass.make(op, [], pref);
n = length(f);
if ( rem(n,2) )
    n = n+1;
    x = trigpts(n);
    v = op(x);
else
    x = trigpts(n);
    v = get(f, 'values');
end
D = valsDiscretizationclass.diffmat(n, 1);
df = D*v;
df_exact = @(x) - 100*x.*exp(-50*x.^2);
err = df_exact(x) - df;
pass(3) = ( norm(err, inf) < 1e3*vscale(f)*eps );
    

a1 = 4; b1 = 3; a2 = 6; b2 = 4;
op = @(x) cos(a1*pi*sin(b1*pi*x)) + 1i*cos(a2*pi*sin(b2*pi*x));
f = testclass.make(op, [], pref);
n = length(f);
if ( rem(n,2) )
    n = n+1;
    x = trigpts(n);
    v = op(x);
else
    x = trigpts(n);
    v = get(f, 'values');
end
D = valsDiscretizationclass.diffmat(n, 1);
df = D*v;
df_exact = @(x) -pi^2*a1*b1*cos(b1*pi*x).*sin(a1*pi*sin(b1*pi*x)) ...
    - 1i*pi^2*a2*b2*cos(b2*pi*x).*sin(a2*pi*sin(b2*pi*x));
err = df_exact(x) - df;
pass(4) = ( norm(err, inf) < 1e5*vscale(f)*eps );
    

%% 2nd-order derivative (odd n & even n):
op = @(x) exp(cos(4*pi*x))-1;
f = testclass.make(op, [], pref);
n = length(f);
x = trigpts(n);
v = get(f, 'values');
D2 = valsDiscretizationclass.diffmat(n, 2);
df2 = D2*v;
df2_exact = @(x) -16*pi^2*exp(cos(4*pi*x)).*(cos(4*pi*x) + cos(4*pi*x).^2 - 1);
err = df2_exact(x) - df2;
pass(5) = ( norm(err, inf) < 5e5*vscale(f)*eps );
    

n = n+1;
x = trigpts(n);
v = op(x);
D2 = valsDiscretizationclass.diffmat(n, 2);
df2 = D2*v;
err = df2_exact(x) - df2;
pass(6) = ( norm(err, inf) < 5e5*vscale(f)*eps );
    

%% Check higher-order derivatives:

% Odd p:
op = @(x) sin(pi*x);
f = testclass.make(op, [], pref);
n = length(f);
x = trigpts(n);
v = get(f, 'values');
D5 = valsDiscretizationclass.diffmat(n, 5);
df5 = D5*v;
df5_exact = @(x) pi^5*cos(pi*x);
err = df5_exact(x) - df5;
pass(7) = ( norm(err, inf) < 1e3*vscale(f)*eps );

n = n+1;
x = trigpts(n);
v = op(x);
D5 = valsDiscretizationclass.diffmat(n, 5);
df5 = D5*v;
err = df5_exact(x) - df5;
pass(8) = ( norm(err, inf) < 1e3*vscale(f)*eps );
    

% Even p:
op = @(x) sin(pi*x);
f = testclass.make(op, [], pref);
n = length(f);
x = trigpts(n);
v = get(f, 'values');
D6 = valsDiscretizationclass.diffmat(n, 6);
df6 = D6*v;
df6_exact = @(x) -pi^6*sin(pi*x);
err = df6_exact(x) - df6;
pass(9) = ( norm(err, inf) < 1e4*vscale(f)*eps );
    

n = n+1;
x = trigpts(n);
v = op(x);
D6 = valsDiscretizationclass.diffmat(n, 6);
df6 = D6*v;
err = df6_exact(x) - df6;
pass(10) = ( norm(err, inf) < 1e5*vscale(f)*eps );
    

end
