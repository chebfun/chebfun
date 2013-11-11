function fig = displayInfoInit(u, pref)
% Utility routine for displaying iteration progress in the solve functions.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences for what we want to show
damped   = pref.damped;
display  = pref.display;%lower(cheboppref('display'));
plotMode = pref.plotting;%lower(cheboppref('plotting'));

% Do we want to print to the command window?
if ( strcmpi(display,'iter') )
    
    % Show info depending on whether we are running in damped mode or not
    if ( damped )
        initString = ['   Iter.       || update ||      Contra.fact.      ', ...
            'length(update)     stepsize    length(u)'];
    else
        initString = ['   Iter.       || update ||      Contra.fact.      ', ...
            'length(update)    length(u)'];
        
    end
    
    % Print to the command window
    fprintf(initString);
    fprintf('\n-------------------------------------------------------------------------------\n');
    
end

% Do we want to show a plot of the initial guess?
if ~( strcmpi(plotMode,'off') )
    fig = figure('name','BVP solver convergence');
    plot(u,'.-'), title('Initial guess of solution')
else
    % If we did not want to plot, we return an empty output instead of a
    % graphics handle.
    fig = [];

end
