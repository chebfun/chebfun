% Test file for conformal.m.

function pass = test_conformal(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Save the current warning state and then clear it
[lastmsg, lastid] = lastwarn();
lastwarn('');

circle = chebfun('exp(pi*1i*x)','trig');
C = real(circle) + .8i*imag(circle);    

% compare the two algorithms on a simple problem
ctr = .1+.1i;
[f,finv] = conformal(C,ctr);
[f2,f2inv] = conformal(C,ctr,'poly');
z = .2-.1i;
err = abs(f(z)-f2(z));
errinv = abs(finv(z)-f2inv(z));
pass(1) = (err < 1e-4) && (errinv < 1e-4);

% snowflake
C = chebfun('exp(pi*1i*t)*(1+.2*cos(6*pi*t))','trig');
[f,finv] = conformal(C);
Z = .5*exp(1i*pi*(1:100)'/100);
err = norm(Z- finv(f(Z)));
pass(2) = (err < 1e-4);

% tanh function
ff = @(z) atanh(2*z)/1.2;
ffinv = @(w) tanh(1.2*w)/2;
C = ffinv(circle);
[f,finv] = conformal(C);
z = -.1i;
pass(3) = abs(f(z)-ff(z)) < 1e-4;
w = .2;
pass(4) = abs(finv(w)-ffinv(w)) < 1e-4;

end
