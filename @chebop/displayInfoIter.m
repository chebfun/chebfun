function displayInfoIter(u, delta, iterNo, normDelta, cFactor, lenDelta, lambda, displayFig, displayTimer, pref)
% Utility routine for displaying iteration progress in the solve functions.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences for what we want to show
damped   = pref.damped;
display  = pref.display;
plotMode = pref.plotting;

% Do we want to print to the command window?
if ( strcmpi(display,'iter') )
    
    % Show info depending on whether we are running in damped mode or not
    if ( damped )
        iterString = sprintf(' %2.2d %13.2e %13.2e    %6i %15.2e ',iterNo, normDelta, ...
            cFactor, lenDelta, lambda);
    else
        iterString = sprintf(' %2.2d %13.2e %13.2e    %6i',iterNo, normDelta, ...
            cFactor, lenDelta);
    end
    
    % Print to the command window
    fprintf(iterString);
    fprintf('\n');
    
end

% Do we want to show a plot of the initial guess?
if ~( strcmpi(plotMode,'off') )
    figure(displayFig);
    subplot(2,1,1)
    plot(u,'.-'), title('Current solution', 'Fontsize', 12)
    subplot(2,1,2)
    plot(delta,'.-'), title('Current correction step', 'Fontsize', 12)
end
