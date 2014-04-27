function displayInfoFinal(u, delta, iterNo, errEstDE, errEstBC, displayFig, displayTimer, pref)
%DISPLAYINFOFINAL
%
% Utility routine for displaying iteration progress in the solve functions. This
% method prints out information after Newton iteration finishes.
%
% See also: displayInfo

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences for what we want to show
display  = pref.display;
plotMode = pref.plotting;

% Do we want to print to the command window?
if ( strcmp(display,'iter') || strcmp(display,'final') )
    % Create a long string of dashes...
    dashString = repmat('-', 1, 62);
    % ... and print to the command window
    fprintf('%s\n', dashString);
    
    if iterNo == 1
        fprintf('Newton''s method converged in 1 iteration\n');
    else
        fprintf('Newton''s method converged in %i iterations.\n',iterNo);
    end
    fprintf(['Final error estimate: %.2e (differential equation) \n' ...
        '%30.2e (boundary conditions). \n\n'], errEstDE, errEstBC);
end

% Do we want to show a plot of the final solution and correction step?
if ~( strcmpi(plotMode,'off') )
    figure(displayFig);
    subplot(2, 1, 1)
    plot(chebfun(u), '.-')
    title('Solution at the end of Newton iteration')
    subplot(2,1,2)
    plot(chebfun(delta),'.-')
    title('Final correction step')
    
end

