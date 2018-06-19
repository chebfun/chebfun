function pass = test_nodots(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

%% Scalar-valued functions

x = chebfun('x');
tol = 1e-14;
    
try 
    a = x*x;
    b = chebfun(@(x) x*x);
    c = chebfun('x*x');
    pass(1) = norm(a - b) + norm(a - c) < tol;
catch
    pass(1) = false;
end

try 
    a = x/(2+x);
    b = chebfun(@(x) x/(2+x));
    c = chebfun('x/(2+x)');
    pass(2) = norm(a - b) + norm(a - c) < tol;
catch
    pass(2) = false;
end

try 
    a = x^2;
    b = chebfun(@(x) x^2);
    c = chebfun('x^2');
    pass(3) = norm(a - b) + norm(a - c) < tol;
catch
    pass(3) = false;
end

%% Vector-valued functions

x = chebfun('x');
x = [x x];
tol = 1e-14;
    
try 
    a = x*x;
    b = x.*x;
    pass(4) = norm(a - b) < tol;
catch
    pass(4) = false;
end

try 
    a = x/(2+x);
    b = x./(2+x);
    pass(5) = norm(a - b) < tol;
catch
    pass(5) = false;
end

try 
    a = x^2;
    b = x.^2;
    pass(6) = norm(a - b) < tol;
catch
    pass(6) = false;
end

%% Test vectorcheck 'off' option

try
    chebfun('x*x', 'vectorcheck', 'off')
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'MATLAB:innerdim');
end

end
