function g = sqrt(f, pref)
% [TODO]: Merge from chebfun-sqrt branch

if ( nargin == 1 )
    pref = chebtech.pref();
end

g = compose(f, @sqrt, [], pref);

end