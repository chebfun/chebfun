function pass = test_galleryballfun( ) 
 
names = {'solharm'}; 
 
N = length(names); 
% Test construction of each gallery function. 
for k = 1:N 
    pass(k) = doesNotCrash(names{k}); 
end 
 
 
if (nargout > 0) 
    pass = all(pass(:)); 
end 
end 
 
function pass = doesNotCrash(name) 
try 
    fn = cheb.galleryballfun(name);  % Test returning the function 
    pass = true; 
catch ME %#ok<NASGU> 
    pass = false; 
end 
end