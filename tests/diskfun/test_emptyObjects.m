function pass = test_emptyObjects( pref ) 
% For empty diskfun objects, does each command deal with them
% appropriately?

% Check things work for empty diskfuns.
f = diskfun();
try
    f + f;
    2*f;
    f*2;
    f.^2;
    2.^f;
    f.^f;
    sqrt(f);
    sum(f);
    integral2(f);
    norm(f);
    squeeze(f);
    diff(f);
    cos(f);
    sin(f);
    sinh(f);
    f.^f + f;
    diag(f);
    trace(f);
    mean(f);
    minandmax2(f);
    median(f);
    flipud(f);
    flipdim(f,1);
    pass = 1;
catch
    pass = 0;
end
end
