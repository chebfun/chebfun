% Test file for @deltafun/isequal.m

function pass = test_isequal(pref)

if (nargin < 1)
    pref = chebfunpref();
end
%%
d1 = deltafun();
d2 = deltafun();

pass(1) = isequal(d1, []) && isequal( [], d1) && isequal( d1, d2); 

f = fun.constructor(@(x) exp(x));
mag = rand(5,5);
loc = rand(1,5);

d1 = deltafun(f, struct('deltaMag', mag, 'deltaLoc', loc));
d2 = deltafun(f, struct('deltaMag', mag, 'deltaLoc', loc));

pass(2) = isequal(d1, d2);

pass(3) = ~isequal( d1, .992312341234*d2 );

d2.deltaMag = d2.deltaMag(1:end-1, :);
pass(4) = ~isequal( d1, d2);

d1.deltaMag = d1.deltaMag(:, 1:end-1);
d1.deltaLoc = d1.deltaLoc(1:end-1);

pass(5) = ~isequal(d1, d2);

f = bndfun(@sin);
d1 = deltafun(f, struct('deltaMag', 1, 'deltaLoc', 0));
d2 = deltafun(f, struct('deltaMag', [1; 0], 'deltaLoc', 0));
pass(6) = isequal(d1, d2);

d1 = deltafun(f, []);
d2 = deltafun(bndfun([]), []);

pass(7) = ~isequal(d1, d2);

end
