function pass = test_galleryballfun( ) 
 
names = {'random', 'zero'}; 
 
N = length(names); 
S = [10,11,12]; 
% Test construction of each gallery function. 
for k = 1:N 
    pass(k) = doesNotCrash(names{k}, S); 
end 
 
 
if (nargout > 0) 
    pass = all(pass(:)); 
end 
end 
 
function pass = doesNotCrash(name, S) 
try 
    fn = cheb.galleryballfun(name, S);  % Test returning the function 
    pass = true; 
catch ME %#ok<NASGU> 
    pass = false; 
end 
end