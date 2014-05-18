function pass = test_chebfun2_plotting( pref )
% Check that the very basic plotting commands do not crash. 
% Alex Townsend, March 2013. 

f = chebfun2(@(x,y) exp(cos(10*x.*y))); 
try 
    plot(f), j = ishold;  
    surf(f), j=j+ishold;
    contour(f),j=j+ishold;
    waterfall(f),j=j+ishold;
    surf(f), j=j+ishold; 
    close all 
    if j==0
        pass(1)=1;
    else
        pass(1) = 0; 
    end
catch
    close all 
    pass(1) = 0;
end

% rank-1 and off [-1 1 -1 1]. 
f = chebfun2(@(x,y) x.*y,[-1 2 -1 2]); 
try 
    plot(f) 
    surf(f)
    contour(f)  
    waterfall(f)
    close all 
    pass(2) = 1; 
catch
    close all 
    pass(2) = 0;
end

try
    % This test comes from PrettyFunctions.m 
h = @(x,y) .75*exp(-((9*x-2).^2 + (9*y-2).^2)/4) + ...
           .75*exp(-((9*x+1).^2)/49 - (9*y+1)/10) + ...
           .5*exp(-((9*x-7).^2 + (9*y-3).^2)/4) - ...
           .2*exp(-(9*x-4).^2 - (9*y-7).^2);
       H = chebfun2(h);
       waterfall(H,'-')
       close all
       pass(3) = 1; 
catch 
    close all
    pass(3) = 0; 
end
end



