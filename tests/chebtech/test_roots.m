% Test file for chebtech/roots.m

function pass = test_roots(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %% Test roots of a bessel function:
    map = @(x) (x+1)*50;
    f = testclass.make(@(x) besselj(0, map(x)), [], [], pref);
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
    pass(n, 1) = norm(r-exact,Inf) < length(f)*f.epslevel;

    %% Test roots of an oscillatory function:
    k = 500;
    f = testclass.make(@(x) sin(pi*k*x), [], [], pref);
    r = roots(f);
    pass(n, 2) = norm(r-(-k:k)'/k, inf) < length(f)*f.epslevel;

    %% Test a perturbed polynomial:
    f = testclass.make( @(x) (x-.1).*(x+.9).*x.*(x-.9) + 1e-14*x.^5, ...
        [], [], pref);
    r = roots(f);
    pass(n, 3) = length(r) == 4 && norm(feval(f, r), inf) < ...
        10*length(f)*f.epslevel;

    %% Test a some simple polynomials:
    f = testclass.make([-1 ; 1], [], [], pref);
    r = roots(f);
    pass(n, 4) = all( r == 0 );

    f = testclass.make([1 ; 0 ; 1]);
    r = roots(f);
    pass(n, 5) = numel(r) == 2 && (norm(r, inf) < f.epslevel);

    %% Test some complex roots:
    f = testclass.make(@(x) 1 + 25*x.^2, [], [], pref);
    r = roots(f, 'complex', 1);

    pass(n, 6) = norm( r - [1i ; -1i]/5, inf) < f.epslevel;

    f = testclass.make(@(x) (1 + 25*x.^2).*exp(x), [], [], pref);
    r = roots(f, 'complex', 1, 'prune', 1);

    pass(n, 7) = norm( r - [1i ; -1i]/5, inf) < 10*length(f)*f.epslevel;

    f = testclass.make(@(x) sin(100*pi*x));
    r1 = roots(f, 'complex', 1, 'recurse', 0);
    r2 = roots(f, 'complex', 1);
    pass(n, 8) = numel(r1) == 383 & numel(r2) >= 383;

    %% Test an array-valued function:
    f = testclass.make(@(x) [sin(pi*x), cos(pi*x)], [], [], pref);
    r = roots(f);
    r2 = [-1 0 1 -.5 .5 NaN].';
    pass(n, 9) = all( r(:) - r2 < 10*length(f)*max(f.epslevel) | isnan(r2) );
end

end
