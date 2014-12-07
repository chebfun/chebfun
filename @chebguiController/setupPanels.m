function handles = setupPanels(handles)
%SETUPPANELS    Populate the panels on CHEBGUI

% Call the subfunctions to populate the panels.
handles = chebguiController.setupPanelInput(handles);
handles = chebguiController.setupPanelType(handles);
handles = chebguiController.setupPanelDiscretization(handles);
handles = chebguiController.setupPanelIVPsolver(handles);

end