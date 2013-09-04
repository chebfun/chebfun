function f = sqrt(f, pref)
% [TODO]: Merge from chebfun-sqrt branch

if ( nargin == 1 )
    pref = fun.pref();
end

pref = f.onefun.pref(pref, pref.fun);

f.onefun = sqrt(f.onefun);

end