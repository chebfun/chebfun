function pass = test_merge(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

d = domain.merge([-inf 1 10 inf], [-inf 1 5 inf]);
pass = all(d == [-Inf     1     5    10   Inf]);


end