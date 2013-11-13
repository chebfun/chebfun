function displayInfoFinal(u, delta, iterNo, errEstDE, errEstBC, displayFig, displayTimer, pref)
% Utility routine for displaying iteration progress in the solve functions.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences for what we want to show
display  = pref.display;
plotMode = pref.plotting;

% Do we want to print to the command window?
% if ( strcmpi(display, 'iter') || strcmpi(display, ' final') )
%     
%     % Show info depending on whether we are running in damped mode or not
%     if ( damped )
%         iterString = sprintf(' %2.2d %12.2e %12.2e    %6i %15.2e %7i ',iterNo, normDelta, ...
%             cFactor, lenDelta, lambda,lenu);
%     else
%         iterString = sprintf(' %2.2d %13.2e %13.2e    %6i',iterNo, normDelta, ...
%             cFactor, lenDelta);
%     end
%     
%     % Print to the command window
%     fprintf(iterString);
%     fprintf('\n');
%     
% end
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

