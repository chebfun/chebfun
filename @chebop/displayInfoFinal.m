function displayInfoFinal(u, delta, iterNo, errEstDE, errEstBC, displayFig, displayTimer, pref)
% Utility routine for displaying iteration progress in the solve functions.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences for what we want to show
display  = pref.display;
plotMode = pref.plotting;

% Do we want to print to the command window?
if strcmp(display,'iter') || strcmp(display,'final')
    fprintf('-------------------------------------------------------------------------------\n');
    if iterNo == 1
        fprintf('Newton''s method converged in 1 iteration\n');
    else
        fprintf('Newton''s method converged in %i iterations.\n',iterNo);
    end
    fprintf('Final error estimate: %.2e (differential equation) \n %29.2e (boundary conditions). \n\n',errEstDE, errEstBC);
end


% Do we want to show a plot of the final solution and correction step?
if ~( strcmpi(plotMode,'off') )
    figure(displayFig);
    subplot(2,1,1)
    plot(u,'.-'), title('Final solution')
    subplot(2,1,2)
    plot(delta,'.-'), title('Final correction step')
    
end

