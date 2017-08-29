function handles = drawLogo(handles)
%DRAWLOGO    Create an axes object in CHEBGUI, and draw the Chebfun logo

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Create a new axes object, store in HANDLES. Hide it until we've finished
% plotting:
handles.fig_logo = axes('Parent', handles.panel_buttons, ...
    'Units', 'normalized', 'Position', [0.025 0.875 0.95 0.1],...
    'Visible', 'off') ;

% Switch focus to our new axes:
axes(handles.fig_logo)

% Load the logo and draw:
logoMat = imread(fullfile(chebfunroot(),'chebguiDemos','chebfunLogo.png'));
image(logoMat)

% Make sure that the logo is displayed:
set(handles.fig_logo, 'Visible', 'on');
axis off
% The following fixes the aspect ratio of the Chebfun logo, even as the CHEBGUI
% window gets rescaled.
daspect([1 1 1])

end
