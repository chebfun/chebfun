function pass = test_surf() 
% Test surf

f = chebfun3(@(x,y,z) cos(x.*y.*z));

pass = 1; 
try 
   surf(f)
catch ME 
    pass = 0; 
end
close all

end