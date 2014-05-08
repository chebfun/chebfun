function pass = test_permute(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

f = chebfun(@sin, pref);

pass(1) = norm(f - permute(f, [1 2])) == 0;

pass(2) = norm(f.' - permute(f, [2 1])) == 0;

try
    permute(f, [1 1]);
    pass(3) = false;
catch 
    pass(3) = true;
end

try
    permute(f, [1 3]);
    pass(4) = false;
catch 
    pass(4) = true;
end

end
    


