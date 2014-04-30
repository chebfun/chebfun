function displayInfoLinear(u, normRes, pref)
%DISPLAYINFOLINEAR
%
% Utility routine for displaying iteration progress in the solve functions. This
% method prints out information for linear problems that are solved using the
% CHEBOP framework.
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
    
    fprintf('Linear equation detected. Converged in one step.\n');
    fprintf('Length of solution: %i.\n',length(chebfun(u)));
    fprintf('Norm of residual: %.2e.\n', normRes);
    
    % Repeat the long string of dashes:
    fprintf('%s\n', dashString);

end

% Do we want to show a plot of the solution?
if ~( strcmpi(plotMode,'off') )
    figure
    plot(chebfun(u), '.-')
    title('Solution of linear ODE')    
end

