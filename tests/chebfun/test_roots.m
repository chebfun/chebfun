function pass = test_roots(pref)

if ( nargin == 0 )
    pref = chebpref();
end

% Test a Sine function:
M = 1000;
f = chebfun(@(x) sin(M*pi*x), [0 1], pref);
exact = linspace(0, 1, M+1).';
r = roots(f);
pass(1) = length(r) == M+1 && norm(exact-r, inf) < epslevel(f);

% Test a polynomial:
p = chebfun( '(x-.1).*(x+.9).*x.*(x-.9) + 1e-14*x.^5' );
exact = roots(p);
pass(2) = length(exact) == 4 && norm(feval(p,exact),inf) < 100*epslevel(p);

% No roots should be returned:
Fs = chebfun({-1, 2}, [-2, 0, 1]);
Fh = chebfun({-1, 0, 1}, [-2, -1, 0, 2]);
rs = roots(Fs, 'nojump');
rh = roots(Fh, 'nojump', 'nozerofun');
pass(3) = isempty(rs) & isempty(rh);

% 0 should be returned as a root in each case:
rs = roots(Fs);
rh = roots(Fh, 'nozerofun');
pass(4) = rs == 0 && numel(rh) == 2 && rh(1) == -1 && rh(2) == 0;

% Test an array-valued function:
f = chebfun(@(x) [sin(2*pi*x), sign(x), x.^2-.5, 1+0*x], [-1, 0, 1], 'extrapolate', 'on');
f.impulses(3,4) = 0;
exact = NaN(5,4);
exact(:,1) = linspace(-1,1,5); exact(1,[2,4]) = [0,1]; exact([1,2],3) = [-1,1]./sqrt(2);
r = roots(f);
pass(5) = all(size(r) == [5,4]) && max(abs(exact(:)-r(:))) < epslevel(f);

% Test roots of a Bessel function on [0, 100]:
f = chebfun(@(x) besselj(0,x),[0 100]);
r = roots(f);
exact = [ 2.40482555769577276862163; 5.52007811028631064959660
8.65372791291101221695437; 11.7915344390142816137431
14.9309177084877859477626; 18.0710639679109225431479
21.2116366298792589590784; 24.3524715307493027370579
27.4934791320402547958773; 30.6346064684319751175496
33.7758202135735686842385; 36.9170983536640439797695
40.0584257646282392947993; 43.1997917131767303575241
46.3411883716618140186858; 49.4826098973978171736028
52.6240518411149960292513; 55.7655107550199793116835
58.9069839260809421328344; 62.0484691902271698828525
65.1899648002068604406360; 68.3314693298567982709923
71.4729816035937328250631; 74.6145006437018378838205
77.7560256303880550377394; 80.8975558711376278637723
84.0390907769381901578795; 87.1806298436411536512617
90.3221726372104800557177; 93.4637187819447741711905
96.6052679509962687781216; 99.7468198586805964702799
];
pass(6) = length(r) == 32 && norm(exact - r, inf) < hscale(f)*epslevel(f);

% Test roots of a complex array-valued CHEBFUN:
f = chebfun(@(x) [1i*x, x, sin(2*pi*x)]);
exact = [[0,0;NaN(4,2)], linspace(-1,1,5).'];
r = roots(f);
pass(7) = all(size(exact) == [5,3]) && max(abs(exact(:) - r(:))) < 10*epslevel(f);

%% Test on singular function: piecewise smooth chebfun - splitting on.

% Set the domain
dom = [-2 7];

pow = -0.5;
op = @(x) (x-dom(1)).^pow.*cos(30*x);
pref.singPrefs.exponents = [pow 0];
pref.enableBreakpointDetection = 1;
f = chebfun(op, dom, pref);
r = roots(f);
r_exact = (((-19:66)+1/2)*pi/30).';
err = r - r_exact;
pass(8) = (norm(err, inf) < 5*get(f, 'vscale')*get(f, 'epslevel'));
    
end