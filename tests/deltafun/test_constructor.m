% Test file for @deltafun/deltafun.m (constructor)

function pass = test_constructor(pref)

if ( nargin < 1 )
     pref = chebpref();
end

%%
deltaTol = pref.deltaPrefs.deltaTol;
proximityTol = pref.deltaPrefs.proximityTol;

d = deltafun();
pass(1) = isempty(d) && isa(d, 'deltafun');

d = deltafun(fun.constructor(0));
pass(2) = isa(d, 'deltafun') && isempty(d.deltaMag) ...
    && isempty(d.deltaLoc) && iszero(d.funPart);

d = deltafun(fun.constructor(@(x) sin(x), [0, 1]));
pass(3) = isa(d, 'deltafun') && isempty(d.deltaMag) ...
    && isempty(d.deltaLoc) && ~iszero(d.funPart);

d = deltafun(rand(3,3), rand(1,3) );
pass(4) = isempty(d.funPart) && all(size(d.deltaMag) == [3, 3]) ...
    && all(size(d.deltaLoc) == [1, 3]);

d = deltafun(rand(5,5), rand(5,1));
pass(5) = isempty(d.funPart) && all(size(d.deltaMag) == [5, 5]) ...
    && all(size(d.deltaLoc) == [1, 5]);

deltas = rand(5,5);
locs = linspace(-.9, .9, 5);
a = randn; b = randn;
f = fun.constructor(@(x) 20*a*sin(10*b*x), [-1,1] );
d = deltafun(f, deltas, locs ); 
pass(6) = iszero(f - d.funPart) && norm(d.deltaMag - deltas, inf) == 0 && ...
    norm(d.deltaLoc(:) - locs(:) ) == 0; 

d = deltafun(f, deltas, locs, pref); 
pass(7) = iszero(f - d.funPart) && norm(d.deltaMag - deltas, inf) == 0 && ...
    norm(d.deltaLoc(:) - locs(:) ) == 0; 


d = deltafun(f, [-1 deltaTol/2 1 deltaTol/2], [-1, 0, .5, 1] );
pass(8) = (norm(d.deltaMag - [-1, 1], inf) == 0) && ...
    (norm(d.deltaLoc - [-1, .5], inf) == 0);

end