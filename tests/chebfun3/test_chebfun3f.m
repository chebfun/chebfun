function pass = test_chebfun3f(pref)
% Test the alternative Chebfun3f constructor for Chebfun3.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

f = @(x,y,z) cos(2*pi*x.*y.*z);
tic()
cf3 = chebfun3(f)
tcf3 = toc()
tic()
cf3f = chebfun3(f,'chebfun3f')
tcf3f = toc()


for samples = 1:1000
    x = rand(1)*2-1;
    y = rand(1)*2-1;
    z = rand(1)*2-1;
    valf = f(x,y,z);
    valcf3 = cf3(x,y,z);
    valcf3f = cf3f(x,y,z);
    err(samples) = abs(valf-valcf3);
    errf(samples) = abs(valf-valcf3f);
end

errChebfun3 = max(err)
errChebfun3f = max(errf)

end