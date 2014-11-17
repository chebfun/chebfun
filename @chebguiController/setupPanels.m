function handles = setupPanels(handles)
%SETUPPANELS    Populate the panels on CHEBGUI

handles = chebguiController.setupPanelInput(handles);
handles = chebguiController.setupPanelType(handles);
handles = chebguiController.setupPanelIVPsolver(handles);
end