function varargout = displayBVPinfo(handles, mode, varargin)

switch mode
    case 'init'
        [displayFig, displayTimer] = displayBVPinfoInit(handles, varargin{:});
        varargout{1} = displayFig;
        varargout{2} = displayTimer;
    case 'iter'
        varargout{1} = displayBVPinfoIter(handles, varargin{:});
    case 'final'
        displayBVPInfoFinal(handles, varargin{:});
    case 'linear'
        displayBVPinfoLinear(handles, varargin{:});
end

end

function [displayFig, displayTimer] = displayBVPinfoInit(handles, u, pref)

% Update the iteration information header on the GUI
initString = ' Iter.       || du ||        Contra.fact.     stepsize        len(du)          len(u)';
set(handles.iter_text,'String', initString);

% Clear the update plot
axes(handles.fig_norm)
plot(0,0)
xlim(handles.xLim)
cla

% Switch focus to the fig_sol plot on the GUI
axes(handles.fig_sol)
% Plot initial guess
plot(chebfun(u), '.-')
title('Initial guess of solution')

% Do different things for the axes depending on whether the solution is real or
% not.
% TODO: This should be calling isreal of chebmatrices.
if isreal(u.blocks{1})
    xlim(handles.xLim)
else
    axis equal
end

drawnow

displayFig = 4;
displayTimer = 5;
end


function displayTimer = displayBVPinfoIter(handles, u, delta, iterNo, normDelta, cFactor, lenDelta, lambda, lenu, displayFig, displayTimer, pref)

% Create a string for displaying information about the iteration.
if ( lambda == 1 )
    iterString = sprintf(' %2.2d %10.2e %10.2e %9.4f %7i %9i ', ...
        iterNo, normDelta, cFactor, lambda, lenDelta, lenu);
else
    iterString = sprintf(' %2.2d %10.2e %10.2e %10.2e %6i %9i ', ...
        iterNo, normDelta, cFactor, lambda, lenDelta, lenu);
end

% Update the iteration information on the GUI
currString = get(handles.iter_list,'String');
set(handles.iter_list,'String', [currString;iterString]);
set(handles.iter_list,'Value',iterNo);


% Switch focus to the fig_sol plot on the GUI
axes(handles.fig_sol)
% Plot initial guess
plot(chebfun(u) ,'.-')
title('Current solution')
% Do different things for the axes depending on whether the solution is real or
% not.
% TODO: This should be calling isreal of chebmatrices.
if isreal(u.blocks{1})
    xlim(handles.xLim)
else
    axis equal
end

% Now plot the Newton updates
% Switch focus to the fig_norm plot on the GUI
axes(handles.fig_norm)
% Plot initial guess
plot(chebfun(delta),'.-')
title('Current correction step')
% Do different things for the axes depending on whether the solution is real or
% not.
% TODO: This should be calling isreal of chebmatrices.
if isreal(u.blocks{1})
    xlim(handles.xLim)
else
    axis equal
end


drawnow


% If the user has pressed the pause button on the GUI, we pause
if ~isempty(handles) && strcmpi(get(handles.button_clear,'String'),'Continue')
    waitfor(handles.button_clear,'String')
end

end

function displayBVPInfoFinal(handles, u, delta, iterNo, errEstDE, errEstBC, displayFig, ...
    displayTimer, pref) %#ok<INUSL>
%DISPLAYINFOFINAL   Utility routine for displaying nonlinear solve progress.
%  This method prints out information after Newton iteration finishes when
%  problems are solved using CHEBGUI.
%
% See also: displayInfo

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% How many steps did we need
if ( iterNo == 1 )
    finalStr = {'Newton''s method converged in 1 iteration\n'};
else
    finalStr = {sprintf('Newton''s method converged in %i iterations.\n',....
        iterNo)};
end

% Show what discretization was used
if ( strcmpi(func2str(pref.discretization), 'ultraS') )
    discString = 'Ultraspherical';
else
    discString = 'Collocation';
end
finalStr = [finalStr; sprintf('Discretization method used: %s. \n', ...
    discString)];
    
% Print info about the final error estimates.
finalStr = [finalStr; ...
    sprintf('Final error estimate: %.2e (differential equation) \n', errEstDE)];
% and the error in the boundary conditions.
finalStr = [finalStr; ...
    sprintf('%30.2e (boundary conditions).', errEstBC)];

% Update the iteration information on the GUI
currString = get(handles.iter_list,'String');
set(handles.iter_list,'String', [currString; finalStr]);

% Set focus to bottom of printed list
set(handles.iter_list,'Value',iterNo + 4);

end

function displayBVPinfoLinear(handles, u, nrmRes, pref)

str = {'Linear equation detected. Converged in one step.'};

% Show what discretization was used
if ( strcmpi(func2str(pref.discretization), 'ultraS') )
    discString = 'Ultraspherical';
else
    discString = 'Collocation';
end
str = [str; sprintf('Discretization method used: %s. \n', discString)];

str = [str; sprintf('Length of solution: %i.\n',length(chebfun(u)))];
str = [str; sprintf('Norm of residual: %.2e.\n', nrmRes)];



set(handles.iter_list,'String', str);
set(handles.iter_list,'Value',1);

end