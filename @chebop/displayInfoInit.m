function [displayFig, displayTimer] = displayInfoInit(u, pref)
% Utility routine for displaying iteration progress in the solve functions.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences for what we want to show
display  = pref.display;
plotMode = pref.plotting;

% Do we want to print to the command window?
if ( strcmpi(display, 'iter') )
    
    % Show info depending on whether we are running in damped mode or not
    initString = ['Iter.   || du ||   Contra.fact.   ', ...
            'stepsize   len(du)   len(u)'];

    % Print to the command window
    fprintf(initString);

    % Create a long string of dashes...
    dashString = repmat('-',1,62);
    % ... and print to the command window
    fprintf('\n%s\n', dashString);
    
end

% Do we want to show a plot of the initial guess?
if ~( strcmpi(plotMode, 'off') )
    displayFig = figure('name','BVP solver convergence');
    plot(u,'.-'), 
    title('Initial guess of solution', 'Fontsize', 12)
else
    % If we did not want to plot, we return an empty output instead of a
    % graphics handle.
    displayFig = [];
end

% Start a timer, which is used to control the pause between successive pauses
displayTimer = tic;
