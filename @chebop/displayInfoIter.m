function displayInfoIter(u, iterNo, normDelta, lenDelta, displayFig, displayTimer, pref)
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
        iterString = sprintf(' %2.2d %13.2e   %25i',iterNo, normDelta, lenDelta);
    else
        iterString = sprintf(' %2.2d %13.2e   %25i',iterNo, normDelta, lenDelta);
    end
    
    % Print to the command window
    fprintf(iterString);
    fprintf('\n');
    
end

% Do we want to show a plot of the initial guess?
if ~( strcmpi(plotMode,'off') )
    figure(displayFig);
    plot(u,'.-'), title('Current solution')
end
