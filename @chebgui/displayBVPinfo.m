function varargout = displayBVPinfo(handles, mode, varargin)
%DISPLAYBVPINFO   Show information on the CHEBGUI figure when solving BVPs.
%
% Calling sequence:
%   VARARGOUT = DISPLAYBVPINFO(HANDLES, MODE, VARARGIN)
% where
%   VARARGOUT when MODE = 'INIT':
%       VARARGOUT{1}:   Dummy output, required for consistency with the
%                       chebop/displayInfoInit method.
%       VARARGOUT{2}:   A tic/toc timer.
%   VARARGOUT when MODE = 'ITER':
%       VARARGOUT{1}:   A tic/toc timer.
%       VARARGOUT{2}:   Boolean variable, equal to 1 if the user pressed the
%                       'Stop' button on the GUI. 
%   HANDLES:            The MATLAB handle to the CHEBGUI figure.
%   VARARGIN:           Useful input arguments for showing information, further
%                       described in 'help chebop/displayInfo'.
%
% See also: chebop/displayInfo

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Go into different modes, depending on where in the Newton iteration we are.
switch mode
    case 'exactInitial'
        displayBVPinfoExactInitial(handles, varargin{:});
    case 'init'
        [displayFig, displayTimer] = displayBVPinfoInit(handles, varargin{:});
        varargout{1} = displayFig;
        varargout{2} = displayTimer;
    case 'iter'
        [varargout{1}, varargout{2}] = displayBVPinfoIter(handles, varargin{:});
    case 'final'
        displayBVPInfoFinal(handles, varargin{:});
    case 'linear'
        displayBVPinfoLinear(handles, varargin{:});
end

end

function displayBVPinfoExactInitial(handles, pref)
%DISPLAYINFOEXACTINITIAL    Called if we pass CHEBOP/SOLVEBVP an initial guess
%                           that solves the BVP.

% Show that we got passed an initial guess that solves the BVP:
str = {'Initial guess appears to be a solution of the BVP.';
    'Returning initial guess as the solution.'};

% Update the information on the GUI:
set(handles.iter_list, 'String',  str);
set(handles.iter_list, 'Value', 1);

end


function [displayFig, displayTimer] = displayBVPinfoInit(handles, u, pref)
%DISPLAYBVPINFOINIT   Show information at the start of Newton iteration
%
% See also: chebop/displayInfoInit

% Update the iteration information header on the GUI
initString = [' Iter.       || du ||        Contraction      stepsize     ',...
    '   len(du)          len(u)'];
set(handles.iter_text, 'String', initString);

% Clear the update plot
axes(handles.fig_norm)
plot(0, 0)
xlim(handles.xLim)
cla

% Switch focus to the fig_sol plot on the GUI
axes(handles.fig_sol)

% Plot initial guess
plot(chebfun(u), '.-', 'linewidth', 2)
set(handles.panel_figSol, 'title', 'Initial guess of solution')

% Remove title from bottom plot panel
set(handles.panel_figNorm, 'title', 'Current correction step')

% Update the fontsize of plots
set(handles.fig_sol, 'fontsize', handles.fontsizePanels);
set(handles.fig_norm, 'fontsize', handles.fontsizePanels);

% Do different things for the axes depending on if the solution is real.
if ( isreal(chebfun(u)) )
    xlim(handles.xLim)
else
    axis equal
end

drawnow

% This only needs to be a dummy output, since the HANDLES variable contains
% pointers to the plots we want to draw on:
displayFig = 42;

% Start a timer, used to control pausing between plotting:
displayTimer = tic;

end

function [displayTimer, stopReq] = displayBVPinfoIter(handles, u, delta, ...
    iterNo, normDelta, cFactor, lenDelta, lambda, lenu, displayFig, ...
    displayTimer, pref) %#ok<INUSD>
%DISPLAYBVPINFOITER   Show information and plot during Newton iteration.
%
% See also: chebop/displayInfoIter

% Create a string for displaying information about the iteration.
if ( lambda == 1 )
    iterString = sprintf(' %2.2d %10.2e %10.2e %9.4f %7i %9i ', ...
        iterNo, normDelta, cFactor, lambda, lenDelta, lenu);
else
    iterString = sprintf(' %2.2d %10.2e %10.2e %10.2e %6i %9i ', ...
        iterNo, normDelta, cFactor, lambda, lenDelta, lenu);
end

% Update the iteration information on the GUI.
currString = get(handles.iter_list, 'String');
set(handles.iter_list, 'String', [currString ; iterString]);
set(handles.iter_list, 'Value', iterNo);

% Switch focus to the fig_sol plot on the GUI.
axes(handles.fig_sol)

% Plot current iterate
plot(chebfun(u), '.-', 'linewidth', 2)
set(handles.panel_figSol, 'title', 'Current solution')

% Do different things for the axes depending on if the solution is real.
if ( isreal(chebfun(u)) )
    xlim(handles.xLim)
else
    axis equal
end

% Now plot the Newton updates.

% Switch focus to the fig_norm plot on the GUI.
axes(handles.fig_norm)

% Plot Newton update
plot(chebfun(delta), '.-',  'linewidth', 2)
set(handles.panel_figNorm, 'title', 'Current correction step')

% Update the fontsize of plots
set(handles.fig_sol, 'fontsize', handles.fontsizePanels);
set(handles.fig_norm, 'fontsize', handles.fontsizePanels);

% Do different things for the axes depending on if the solution is real.
if ( isreal(chebfun(u)) )
    xlim(handles.xLim)
else
    axis equal
end

drawnow

% If the user has pressed the pause button on the GUI, we pause.
if ( ~isempty(handles) ...
            && strcmpi(get(handles.button_clear, 'String'), 'Continue') )
    waitfor(handles.button_clear, 'String')
end

% If the user has pressed the stop button on the GUI, we stop the Newton
% iteration.
if ( ~isempty(handles) ...
            && strcmpi(get(handles.button_solve, 'String'), 'Solve') )
    stopReq = true;
else
    stopReq = false;
end

% Restart the timer:
displayTimer = tic;

end

function displayBVPInfoFinal(handles, u, delta, iterNo, errEstDE, errEstBC, ...
	displayFig, displayTimer, pref) %#ok<INUSL>
%DISPLAYBVPINFOFINAL   Utility routine for displaying nonlinear solve progress.
%
%  This method prints out information after Newton iteration finishes when
%  problems are solved using CHEBGUI.
%
% See also: displayBVPInfo, chebop/displayInfoFinal

% How many steps did we need
if ( iterNo == 1 )
    finalStr = {'Newton''s method converged in 1 iteration\n'};
else
    finalStr = {sprintf('Newton''s method converged in %i iterations.',....
        iterNo)};
end

% Show what discretization was used
if ( strcmpi(func2str(pref.discretization), 'coeffs') || ...
        isequal(pref.discretization, @ultraS) || ...
        isequal(pref.discretization, @trigspec))
    discString = 'Coefficients';
else
    discString = 'Values';
end
finalStr = [finalStr; sprintf('Discretization basis used: %s.', ...
    discString)];
    
% Print info about the final error estimates...
finalStr = [finalStr; ...
    sprintf('Final error estimate: %.2e (differential equation) ', errEstDE)];
% ...and the error in the boundary conditions.
finalStr = [finalStr; ...
    sprintf('%30.2e (boundary conditions).', errEstBC)];

% Update the iteration information on the GUI
currString = get(handles.iter_list, 'String');
set(handles.iter_list,'String', [currString ; finalStr]);

% Set focus to bottom of printed list
set(handles.iter_list, 'Value', iterNo + 4);

end

function displayBVPinfoLinear(handles, u, nrmRes, pref)
%DISPLAYBVPINFOLINEAR   Display methods tailored to linear BVPs.

% Show that we detetected a linear equation:
str = {'Linear equation detected. Converged in one step.'};

% Show what discretization was used:
if ( strcmpi(func2str(pref.discretization), 'coeffs') || ...
        isequal(pref.discretization, @ultraS) || ...
        isequal(pref.discretization, @trigspec))
    discString = 'Coefficients';
else
    discString = 'Values';
end

% Concatenate strings:
str = [str ; sprintf('Discretization basis used: %s.',  discString)];
str = [str ; sprintf('Length of solution: %i.', length(chebfun(u)))];
str = [str ; sprintf('Norm of residual: %.2e.', nrmRes)];

% Update the information on the GUI:
set(handles.iter_list, 'String',  str);
set(handles.iter_list, 'Value', 1);

% Plot
axes(handles.fig_sol)
plot(u, 'linewidth', 2)
% Do different things depending on whether the solution is real or not
if ( isreal(u{1}) )
    axis tight
else
    axis equal
end

set(handles.panel_figSol, 'title', 'Solution')

end
