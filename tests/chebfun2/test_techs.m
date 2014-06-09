function pass = test_techs( prefs )
% Check that Chebfun2 can work with the different techs. 

if ( nargin < 1 )
    pref = chebfunpref(); 
end

try 
    pref.tech = @chebtech2;
    chebfun2(@(x,y) cos(x.*y), pref); 
    pass(1) = 1;
catch
    pass(1) = 0; 
end

try 
    pref.tech = @chebtech1;
    hebfun2(@(x,y) cos(x.*y), pref);
    pass(2) = 1;
catch
    pass(2) = 0; 
end

% try 
%     pref.tech = @fourtech;
%     chebfun2(@(x,y) cos(pi*y).*sin(pi*x), pref);
%     pass(3) = 1;
% catch
%     pass(3) = 0; 
% end

end
