function pass = test_chebfun_plot(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

h = figure();

try 
    x = chebfun(@(x) x, pref);
    S = [-1+1i*x -1i+x 1i*x];  
    plot(S);
    close(h)
    pass(1) = true;
catch ME
    close(h)
    rethrow(ME)
end

end
