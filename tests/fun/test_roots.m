% Test file for fun/roots.m

function pass = test_roots(pref)

if ( nargin < 1 )
    pref = fun.pref;
end

pass = zeros(1, 9); % Pre-allocate pass matrix
for n = 1:1  %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
    else 
        testclass = unbndfun();
    end

    %% Test roots of a bessel function:
    map = @(x) (x+2)*100/9;
    f = testclass.make(@(x) besselj(0, map(x)), dom, [], [], pref);
    r = map(roots(f));
    exact = [   2.40482555769577276862163; 5.52007811028631064959660
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
    pass(n, 1) = norm(r-exact,Inf) < length(f)*get(f, 'epslevel');

    %% Test roots of an oscillatory function:
    k = 100;
%     f = chebtech.constructor(@(x) sin(pi*k*x), [], [], chebtech.pref);
    f = testclass.make(@(x) sin(pi*k*x), dom, [], [], pref);
    r = roots(f);
    
    pass(n, 2) = norm(r-(-2*k:7*k)'/k, inf) < get(f, 'epslevel').*get(f, 'vscale');

    %% Test a perturbed polynomial:
    f = testclass.make( @(x) (x-.1).*(x+.9).*x.*(x-.9) + 1e-14*x.^5, dom, [], [], pref);
    r = roots(f);
    pass(n, 3) = length(r) == 4 && norm(feval(f, r), inf) < 10*get(f, 'epslevel').*get(f, 'vscale');

    %% Test a some simple polynomials:
    f = testclass.make([-2 ; 7], dom, [], [], pref);
    r = roots(f);
    pass(n, 4) = r < get(f, 'epslevel').*get(f, 'vscale');;

    f = testclass.make([20.25 ; 0 ; 20.25]);
    r = roots(f);
    pass(n, 5) = numel(r) == 2 && (norm(r, inf) < get(f, 'epslevel').*get(f, 'vscale'));

    %% Test some complex roots:
    f = testclass.make(@(x) 1 + x.^2, dom, [], [], pref);
    r = roots(f, 'complex', 1);

    pass(n, 6) = norm( r - [1i ; -1i], inf) < get(f, 'epslevel').*get(f, 'vscale');

    f = testclass.make(@(x) (1 + 25*x.^2).*exp(x), [-1 1], [], [], pref);
    r = roots(f, 'complex', 1, 'prune', 1);

    pass(n, 7) = norm( r - [1i ; -1i]/5, inf) < get(f, 'epslevel').*get(f, 'vscale');

    f = testclass.make(@(x) sin(10*pi*x), dom);
    r1 = roots(f, 'complex', 1, 'recurse', 0);
    r2 = roots(f, 'complex', 1);
    pass(n, 8) = numel(r1) == 195 & numel(r2) >= 195;

    %% Test an array-valued function:
    f = testclass.make(@(x) [sin(pi*x), cos(pi*x), x.^2+1], dom, [], [], pref);
    r = roots(f);
    r2 = [-2:7 -1.5:6.5 NaN(1,11)].';
    pass(n, 9) = all( r(:) - r2 < max(get(f, 'epslevel').*get(f, 'vscale')) | isnan(r2) );
end

end