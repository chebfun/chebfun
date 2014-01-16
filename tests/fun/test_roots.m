% Test file for fun/roots.m

function pass = test_roots(pref)

if ( nargin < 1 )
    pref = chebpref();
end

% Set a domain for BNDFUN.
dom = [-2 7];

%% 
% Test roots of a bessel BNDFUN:
map = @(x) (x+2)*100/9;
f = bndfun(@(x) besselj(0, map(x)), dom, [], [], pref);
r = map(roots(f));
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
          96.6052679509962687781216; 99.7468198586805964702799 ];
pass(1) = norm(r-exact,Inf) < length(f)*get(f, 'epslevel');

%% 
% Test roots of an oscillatory BNDFUN:
k = 100;
f = bndfun(@(x) sin(pi*k*x), dom, [], [], pref);
r = roots(f);
pass(2) = norm(r-(-2*k:7*k)'/k, inf) < get(f, 'epslevel').*get(f, 'vscale');

%%
% Test a perturbed polynomial BNDFUN:
f = bndfun( @(x) (x-.1).*(x+.9).*x.*(x-.9) + 1e-14*x.^5, dom, [], [], pref);
r = roots(f);
pass(3) = length(r) == 4 && norm(feval(f, r), inf) < 10*get(f, 'epslevel').*get(f, 'vscale');

%%
%  Test a some simple polynomials BNDFUN:
f = bndfun([-2 ; 7], dom, [], [], pref);
r = roots(f);
pass(4) = r < get(f, 'epslevel').*get(f, 'vscale');;

f = bndfun([20.25 ; 0 ; 20.25]);
r = roots(f);
pass(5) = numel(r) == 2 && (norm(r, inf) < get(f, 'epslevel').*get(f, 'vscale'));

%%
% Test some complex roots of BNDFUN:
f = bndfun(@(x) 1 + x.^2, dom, [], [], pref);
r = roots(f, 'complex', 1);
pass(6) = norm( r - [1i ; -1i], inf) < get(f, 'epslevel').*get(f, 'vscale');

f = bndfun(@(x) (1 + 25*x.^2).*exp(x), [-1 1], [], [], pref);
r = roots(f, 'complex', 1, 'prune', 1);
pass(7) = norm( r - [1i ; -1i]/5, inf) < get(f, 'epslevel').*get(f, 'vscale');

f = bndfun(@(x) sin(10*pi*x), dom);
r1 = roots(f, 'complex', 1, 'recurse', 0);
r2 = roots(f, 'complex', 1);
pass(8) = numel(r1) == 195 & numel(r2) >= 195;

%%
% Test an array-valued BNDFUN:
f = bndfun(@(x) [sin(pi*x), cos(pi*x), x.^2+1], dom, [], [], pref);
r = roots(f);
r2 = [-2:7 -1.5:6.5 NaN(1,11)].';
pass(9) = all( r(:) - r2 < max(get(f, 'epslevel').*get(f, 'vscale')) | isnan(r2) );
    
%% 
% Test on singular BNDFUN.
pow = -0.5;
op = @(x) (x-dom(1)).^pow.*cos(x);
pref.singPrefs.exponents = [pow 0];
f = bndfun(op, dom, [], [], pref);
r = roots(f);
r_exact = [-1/2; 1/2; 3/2]*pi;
err = r - r_exact;
pass(10) = (norm(err, inf) < 1e2*get(f, 'vscale')*get(f, 'epslevel'));
    
%% 
% [TODO]: Run a few tests for UNBNDFUN.
end