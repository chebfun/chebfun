% Test file for @chebfun/legpoly.m.

function pass = test_legpoly(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Check coefficients of Legendre polynomials.
Q = legpoly(4:-1:0);
pass(1) = norm(legpoly(Q) - eye(5), 'fro') < 10*vscale(Q)*epslevel(Q);

% Check coefficients of a linear combination of Legendre polynomials.
Q = legpoly(9:-1:0);
v = (1:1:10).';
Qv = Q*v;
pass(2) = norm(legpoly(Qv) - v.', Inf) < 10*vscale(Qv)*epslevel(Qv);

% Check a less trivial smooth example.
f = chebfun(@(x) [sin(x - 0.3) cos(x - 0.3)]);
c = legpoly(f, 5);
c_exact = [-0.0026889804057628215030252  0.0086927426357443276344859
           -0.0602302090841307493414730 -0.0186313869912884662151106
            0.0916630569532407635738723 -0.2963217439563513424659842
            0.8631522851187122776801783  0.2670042907205986623437124
           -0.2486716793299505096289379  0.8038879363274420091869815].';
pass(3) = norm(c(:) - c_exact(:), Inf) < 10*vscale(f)*epslevel(f);

% Check some piecewise examples.
f = chebfun(@(x) abs(x - 0.3), [-1 0.3 1]);
c = legpoly(f, 5);
c_exact = [-0.0574494375
            0.2173762500
            0.5175625000
           -0.4365000000
            0.5450000000].';
pass(4) = norm(c - c_exact, Inf) < 10*vscale(f)*epslevel(f);

f = chebfun(@(x) [abs(x - 0.3).^3 exp(x)], [-1 -0.5 0.3 0.5 1]);
c = legpoly(f, 5);
c_exact = [ 0.09643353890625000000000  0.009965128148869178524617
           -0.24012341250000000000000  0.070455633668489027815297
            0.78877862500000000000000  0.357814350647372460479053
           -0.71513550000000000000000  1.103638323514326964786571
            0.38702500000000000000000  1.175201193643801456882382].';
pass(5) = norm(c(:) - c_exact(:), Inf) < 10*vscale(f)*epslevel(f);

% Check operation for row CHEBFUNs.
c = legpoly(f.', 5);
pass(6) = norm(c(:) - c_exact(:), Inf) < 10*vscale(f)*epslevel(f);

end
