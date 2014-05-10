function displayInfoFinal(u, delta, iterNo, errEstDE, errEstBC, displayFig, ...
    displayTimer, pref) %#ok<INUSL>
%DISPLAYINFOFINAL   Utility routine for displaying nonlinear solve progress.
%  This method prints out information after Newton iteration finishes and a plot
%  of the final solution obtain, if it is specified in the CHEBOPPREF argument
%  passed to this method that the user wants to see the info/plot.
%
% See also: displayInfo

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Document inputs.

% Obtain preferences for what we want to show
display  = pref.display;
plotMode = pref.plotting;

% Do we want to print to the command window?
if ( strcmp(display,'iter') || strcmp(display,'final') )
    % Create a long string of dashes...
    dashString = repmat('-', 1, 62);
    % ... and print to the command window
    fprintf('%s\n', dashString);
    
    % How many steps did we need
    if ( iterNo == 1 )
        fprintf('Newton''s method converged in 1 iteration\n');
    else
        fprintf('Newton''s method converged in %i iterations.\n', iterNo);
    end
    
    % Print info about the final error estimates.
    fprintf(['Final error estimate: %.2e (differential equation) \n' ...
        '%30.2e (boundary conditions). \n\n'], errEstDE, errEstBC);
end

if ( ~strcmpi(plotMode, 'off') )
    % Plot of the final solution and correction step:
    figure(displayFig);
    subplot(2, 1, 1)
    plot(chebfun(u), '.-')
    title('Solution at the end of Newton iteration')
    subplot(2, 1, 2)
    plot(chebfun(delta), '.-')
    title('Final correction step')
end

end

