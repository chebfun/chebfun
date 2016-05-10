function varargout = displayInfoLinear(u, normRes, pref)
%DISPLAYINFOLINEAR   Utility routine for displaying linear solve progress.
%   DISPLAYINFOLINEAR(U, NORMRES) pretty-prints the CHEBFUN solution U to a
%   linear CHEBOP along with the norm of the residual, NORMRES, to the screen. U
%   is also plotted in a new figure window.
%
%   H = DISPLAYINFOLINEAR(U, NORMRES) returns a figure handle to the plot.
%
%   DISPLAYINFOLINEAR(U, NORMRES, PREF) accepts preferences through a CHEBOPPREF
%   object PREF. In particular, PREF.DISPLAY = 'off' will prevent printing of any
%   information, and PREF.PLOTTING = 'off' prevents plotting.
%
% See also: CHEBOP.DISPLAYINFO.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Obtain preferences for what we want to show:
if ( nargin < 3 )
    pref = cheboppref();
end
display  = pref.display;
plotMode = pref.plotting;

% Do we want to print to the command window?
if ( strcmp(display, 'iter') || strcmp(display, 'final') )
    
    % Create a long string of dashes...
    dashString = repmat('-', 1, 62);
    % ... and print to the command window
    fprintf('%s\n', dashString);

    % Determine the length of the solution:
    if ( isa(u, 'chebmatrix') )
        l = max(cellfun(@length, u.blocks));
    else
        l = length(u);
    end
    
    % Plot info:
    fprintf('Linear equation detected. Converged in one step.\n');
    
    % Show what discretization was used
    if ( strcmpi(func2str(pref.discretization), 'coeffs') || ...
            isequal(pref.discretization, @ultraS) || ...
            isequal(pref.discretization, @trigspec) )
        discString = 'Coefficients';
    else
        discString = 'Values';
    end
    fprintf('Discretization basis used: %s. \n', discString);
    fprintf('Length of solution: %i.\n', l);
    fprintf('Norm of residual: %.2e.\n', normRes);
    
    % Repeat the long string of dashes:
    fprintf('%s\n', dashString);

end

% Do we want to show a plot of the solution?
if ( ~strcmpi(plotMode, 'off') )
    figure
    h = plot(chebfun(u), '.-');
    title('Solution of linear ODE')    
else
    h = [];
end

if ( nargout > 0 )
    varargout = {h};
end

end

