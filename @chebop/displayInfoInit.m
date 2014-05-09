function [displayFig, displayTimer] = displayInfoInit(u, pref)
%DISPLAYINFOINIT   Utility routine for displaying linear solve progress.
% This method prints out information before the Newton iteration starts.
% The outputs of this function are:
%   displayFig:     A handle to the MATLAB figure plots will be drawn on.
%   displayTimer:   A tic/toc timer, so that we can control the pause between
%                   plots.
%
% See also: DISAPLYINFO.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Document inputs.

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
    dashString = repmat('-', 1, 62);
    % ... and print to the command window
    fprintf('\n%s\n', dashString);
    
end

% Do we want to show a plot of the initial guess?
if ( ~strcmpi(plotMode, 'off') )
    displayFig = figure('name', 'BVP solver convergence');
    plot(chebfun(u), '.-'), 
    title('Initial guess of solution', 'Fontsize', 12)
    drawnow
    
    % Do we need to pause?
    if ( strcmpi(plotMode, 'pause') )
        pause
    end
    
else
    % If we did not want to plot, we return an empty output instead of a
    % graphics handle.
    displayFig = [];
    
end

% Start a timer, which is used to control the pause between successive pauses
displayTimer = tic;

end