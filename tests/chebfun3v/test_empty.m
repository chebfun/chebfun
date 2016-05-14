function pass = test_empty()
% Does each Chebfun3 command deal with empty inputs appropriately?

pass = 1; 
F = chebfun3v;
try 
    conj(F);
    cross(F,F);
    F.';
    curl(F);
    divergence(F);
    dot(F,F);
    F(1,1,1);
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