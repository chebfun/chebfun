% Test for trigcf.m.
function pass = test_trigcf(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Set an error tolerance.
tol = 1.0e-12;

%% check for emptiness:
p = trigcf(chebfun);
pass(1) = isempty(p);


%% check for a simple trigfun 
fh = @(x) cos(4*pi*x) + 1;
f = chebfun(fh, [0, 2], 'trig');
[p, s] = trigcf(f, 4);
pass(2) = norm(p-f, inf) < tol && abs(s) < tol;

%%
% Test an example that is not based on the domain [-1, 1].
f = chebfun(@(x) exp(cos(pi*x)), [2 6], 'trig');
p = trigcf(f, 5);
pass(3) = (length(p) == 11) && all(domain(p) == [2, 6]); 


% Test CF with quasimatrix and array-valued input to ensure the output is of
% the correct form.
f = chebfun(@(x) [sin(pi*x) cos(pi*x) exp(sin(pi*x))], [-1 1], 'trig');
[p, s] = trigcf(f, 2);
pass(4) = (numel(p) == 3) && (numel(s) == 3);

f = cheb2quasi(f).';
[p, s] = trigcf(f, 2);
pass(5) = (numel(p) == 3) && (numel(s) == 3) && ...
    p(1,:).isTransposed && iscolumn(s);


f = chebfun(@(x) exp(cos(pi*x) + sin(2*pi*x)), 'trig');
M = ceil((length(f)-1)/2);
pass(6) = norm(trigcf(f,M) - f, inf) < 100*eps*f.vscale;
pass(7) = norm(trigcf(f, M+10) - f, inf) < 100*eps*f.vscale;
f = chebfun(@(x) cos(4*pi*x), 'trig');
[p, s] = trigcf(f, 3);
pass(8) = norm(p, inf) < 100*eps*f.vscale;
pass(9) = abs(s - 1) < 100*eps*f.vscale;

end

