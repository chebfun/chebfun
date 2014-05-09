function displayTimer = displayInfoIter(u, delta, iterNo, normDelta, cFactor, lenDelta, lambda, lenu, displayFig, displayTimer, pref)
%DISPLAYINFOITER
%
% Utility routine for displaying iteration progress in the solve functions. This
% method prints and plots during the Newton iteration.
%
% See also: displayInfo

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences for what we want to show
display  = pref.display;
plotMode = pref.plotting;

% Do we want to print to the command window?
if ( strcmpi(display,'iter') )
    
    % Create a string with information about the convergence. Want a slightly
    % neater output if we are taking a full step.
    if ( lambda == 1 )
        iterString = sprintf(' %2.2d %12.2e %12.2e %11.4f %7i %9i ', ...
            iterNo, normDelta, cFactor, lambda, lenDelta, lenu);        
    else
        iterString = sprintf(' %2.2d %12.2e %12.2e %12.2e %6i %9i ', ...
            iterNo, normDelta, cFactor, lambda, lenDelta, lenu);
    end
    
    % Print to the command window
    fprintf(iterString);
    fprintf('\n');
    
end

% Do we want to show a plot of the initial guess?
if ~( strcmpi(plotMode,'off') )
    
    % Do we need to pause?
    if ( strcmpi(plotMode, 'pause') )
        pause
    elseif ( isnumeric(plotMode) )
        % Specified a numerical value for the pause
        pause( plotMode - toc(displayTimer) )
    end
    
    % Set focus on the correct figure
    figure(displayFig);
    
    % Plot current solution
    subplot(2, 1, 1)
    plot(chebfun(u), '.-')
    title('Current solution', 'Fontsize', 12)
    
    % Plot current correction
    subplot(2, 1, 2)
    plot(chebfun(delta),'.-')
    title('Current correction step', 'Fontsize', 12)
    
    drawnow

    % Restart the timer
    displayTimer = tic;
end
