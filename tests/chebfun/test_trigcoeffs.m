function pass = test_trigcoeffs(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% Test ordering of coefficients is correct
f_test = @(x) 1 + (2+2i)*exp(1i*2*pi*x) + (5+5i)*exp(-1i*5*pi*x);
f = chebfun(f_test, 'trig');
c = trigcoeffs(f);
c_exact = [(5+5i) 0 0 0 0 1 0 (2+2i) 0 0 0].';
err = c-c_exact;
pass(1) = norm(err,inf) < 1e2*vscale(f)*eps;


% Test sin/cos coefficient form is returned correctly.
f_test = @(x) 1 + 5*cos(5*pi*x) + 10*cos(10*pi*x) - 7*sin(7*pi*x) + 8*sin(8*pi*x); 
f = chebfun(f_test, 'trig');
[a,b] = trigcoeffs(f);
a_exact = [1 0 0 0 0 5 0 0 0 0 10].';
b_exact = [0 0 0 0 0 0 -7 8 0 0].';
err = [(a-a_exact);(b-b_exact)];
pass(2) = norm(err,inf) < 1e2*vscale(f)*eps;


% Test on simple combination of Fourier modes on symmetric domain about the
% origin.
f_test = @(x) 2 + 2*cos(2*x) + sin(x);
dom = [-pi pi];
f = chebfun(f_test,dom, 'trig');
c = trigcoeffs(f);
c_exact = [1 0.5i 2 -0.5i 1].';
err = c-c_exact;
pass(3) = norm(err,inf) < 1e2*vscale(f)*eps;


% Now change domains and check that the same result is given for the 
% coefficients.
dom = [0.1 0.1+2*pi];
f = chebfun(f_test,dom, 'trig');
c = trigcoeffs(f);
err = c-c_exact;
pass(4) = norm(err,inf) < 1e2*vscale(f)*eps;


% Check for error when input is a chebfun and not a trigfun.
try
    f = chebfun(@(x) exp(cos(x)),dom);
    c = trigcoeffs(f);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:trigcoeffs:chebfun');
end


end
