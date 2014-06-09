function pass = test_techs( prefs )
% Check that Chebfun2 can work with the different techs. 

if ( nargin < 1 )
    prefs = chebfunpref(); 
end

try 
    chebfunpref.setDefaults('tech',@chebtech2)
    f = chebfun2(@(x,y) cos(x.*y)); 
    pass(1) = 1;
catch
    pass(1) = 0; 
end

try 
    chebfunpref.setDefaults('tech',@chebtech1)
    f = chebfun2(@(x,y) cos(x.*y));
    pass(2) = 1;
catch
    pass(2) = 0; 
end

% try 
%     chebfunpref.setDefaults('tech',@fourtech)
%     f = chebfun2(@(x,y) cos(pi*y).*sin(pi*x));
%     pass(3) = 1;
% catch
%     pass(3) = 0; 
% end

end
