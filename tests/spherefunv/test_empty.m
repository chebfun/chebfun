function pass = test_empty( ) 
% For empty spherfunv objects, does each command deal with them
% appropriately?

%% 
% Check that the spherefunv commands work for empty Spherefunv objects. 

pass = 1; 
F = spherefunv;
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
    vorticity(F);
    normal(F);
catch
   pass = 0;  
end

end
