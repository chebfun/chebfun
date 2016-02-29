% Test file for @trigtech/sample.m

function pass = test_sample(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

f = trigtech(@(x) exp(sin(pi*x - 0.1)));

% Test on a grid equal to length(f).
[v, p] = sample(f);
p_ex = trigtech.trigpts(length(f));
v_ex = feval(f, p_ex);
pass(1) = (norm(p - p_ex) < 100*eps) && (norm(v - v_ex) < 100*eps);

% Test on a grid shorter than length(f).
m = round(length(f)/2);
[v, p] = sample(f, m);
p_ex = trigtech.trigpts(m);
v_ex = feval(f, p_ex);
pass(2) = (norm(p - p_ex) < 100*eps) && (norm(v - v_ex) < 100*eps);

% Test on a grid longer than length(f).
m = round(2*length(f));
[v, p] = sample(f, m);
p_ex = trigtech.trigpts(m);
v_ex = feval(f, p_ex);
pass(3) = (norm(p - p_ex) < 100*eps) && (norm(v - v_ex) < 100*eps);

end

