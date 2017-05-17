function pass = test_empty( ) 
% For empty diskfunv objects, does each command deal with them
% appropriately?

%% 
% Check that the diskfunv commands work for empty diskfunv objects. 

pass = 1; 
F = diskfunv;
try 
    cross(F,F);
    F.';
    curl(F);
    divergence(F);
    dot(F,F);
    F(1,1);
    isempty(F);
    F-F;
    F.*F;
    norm(F);
    F+F;
    F.^2;
    F*F;
    F';
catch
   pass = 0;  
end

end
