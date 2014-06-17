function pass = test_empty( pref ) 
% For empty chebfun2v objects, does each command deal with them
% appropriately?
% Alex Townsend, March 2013. 

%% 
% Check that the chebfun2v commands work for empty Chebfun2v objects. 

pass = 1; 
F = chebfun2v;
try 
    conj(F);
    cross(F,F);
    F.';
    curl(F);
    divergence(F);
    dot(F,F);
    F(1,1);
    imag(F); 
    isempty(F);
    laplacian(F);
    F-F;
    F.*F;
    norm(F);
    F+F;
    F.^2;
    real(F); 
    roots(F);
    F*F;
    F';
catch
   pass = 0;  
end

end
