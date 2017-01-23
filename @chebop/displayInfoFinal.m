function displayInfoFinal(u, delta, iterNo, errEstDE, errEstBC, displayFig, ...
    displayTimer, pref) %#ok<INUSL>
%DISPLAYINFOFINAL   Utility routine for displaying nonlinear solve progress.
%  This method prints out information after Newton iteration finishes and a plot
%  of the final solution obtain, if it is specified in the CHEBOPPREF argument
%  passed to this method that the user wants to see the info/plot.
%
%  Calling sequence:
%       DISPLAYINFOFINAL(U, DELTA, ITERNO, ERRESTDE, ERRESTBC, DISPLAYFIG, ...
%           DISPLAYTIMER, PREF)
%
%  where:
%   U:              Solution of the nonlinear BVP.
%   DELTA:          Final Newton correction
%   ITERNO:         Total number of iterations needed.
%   ERRESTDE:       Estimate of the error of the solution of the differntial eq.
%   ERRESTBC:       Estimate of the error of satisfying the boundary conditions.
%   DISPLAYFIG:     Handle to a Matlab figure for plotting.
%   DISPLAYTIMER:   A tic/toc timer for controlling pauses.
%   PREF:           A CHEBOPPREF() object.
%
% See also: displayInfo

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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
        fprintf('Newton''s method converged in 1 iteration.\n');
    else
        fprintf('Newton''s method converged in %i iterations.\n', iterNo);
    end
    
    % Show what discretization was used
    if ( strcmpi(func2str(pref.discretization), 'coeffs') || ...
            isequal(pref.discretization, @ultraS) || ...
            isequal(pref.discretization, @trigspec) )
        discString = 'Coefficients';
    else
        discString = 'Values';
    end
    fprintf('Discretization basis used: %s. \n', discString);
    
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

