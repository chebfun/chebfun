function pass = test_size() 
% Test SIZE

% Check two components: 
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(x+y+z));
pass(1) = all(size(F) == [2 inf inf inf]);

% Check three components: 
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) cos(z));
pass(2) = all(size(F) == [3 inf inf inf]);

end