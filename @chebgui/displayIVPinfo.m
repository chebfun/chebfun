function displayIVPinfo(handles, u, isIVP, varargin)
%DISPLAYIVPINFO   Show information on the CHEBGUI figure when solving IVPs.
%
% Calling sequence:
%   DISPLAYBVPINFO(HANDLES, VARARGIN)
% where
%   HANDLES:   The MATLAB handle to the CHEBGUI figure.
%   U:         The solution computed.
%   ISIVP:     Equal to 1 if we've just solved an IVP, 0 otherwise (FVP).
%   VARARGIN:  Useful input arguments for showing information, further described
%              in 'help chebop/displayInfo'.
%
% See also: chebop/displayInfo.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Are we dealing with a CHEBFUN or a CHEBMATRIX?
if (isa(u, 'chebfun'))
    len = length(u);
else
    len = max(cellfun(@length, u.blocks));
end

% Did we just solve an IVP or an FVP?
if ( isIVP == 1 )
    detStr = 'Initial value problem detected';
else
    detStr = 'Final value problem detected';
end

% String to be shown in CHEBGUI
str = {detStr, ;
    sprintf('Number of pieces of the solution: %i.', ...
    length(u.domain)-1);
    sprintf('Total length of solution: %i.', len)};
% Show the string:
set(handles.iter_list, 'String',  str)
% Focus on top of the box:
set(handles.iter_list, 'Value', 1)

% Plot
axes(handles.fig_sol)
plot(u, 'Linewidth', 2)
% Do different things depending on whether the solution is real or not
if ( (isa(u, 'chebfun') && isreal(u)) || (isa(u, 'chebmatrix') && isreal(u{1})) )
    xlim(handles.xLim)
else
    axis equal
end

end

