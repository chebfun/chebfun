function pass = test_galleryballfunv( ) 
 
names = {'zero'}; 
 
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
    vn = cheb.galleryballfunv(name);  % Test returning the vector 
    pass = true; 
catch ME %#ok<NASGU> 
    pass = false; 
end 
end