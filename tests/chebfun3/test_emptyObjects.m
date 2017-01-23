function pass = test_emptyObjects() 
% Does each chebfun3 command deal with empty chebfun3 objects, 
% appropriately?

% Check things work for empty chebfun3 objects.
f = chebfun3();
try
    f + f;
    2*f;
    f*2;
    f.^2;
    2.^f;
    f.^f;
    sqrt(f);
    sum(f);
    sum2(f);
    sum3(f);
    integral(f);
    diff(f);
    sin(f);
    cos(f);
    sinh(f);
    f.^f + f;
    tucker(f)
    mean(f);
    max3(f);
    norm(f);    
    minandmax3(f);
    permute(f);
    pass = 1;
catch
    pass = 0;
end
end
