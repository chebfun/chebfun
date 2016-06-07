function handles = changeFontsize(handles, change)

currentFontsize = handles.fontsizePanels;
newFontsize = currentFontsize + change;

% Allow max fontsize 22 and min fontsize 6:
maxFontsize = 22;
minFontsize = 6;

if ( (newFontsize > maxFontsize) || (newFontsize < minFontsize) )
    % Maximum or minimum fontsize reached -- do nothing.
    return
end

% Get the field names in CHEBGUI:
names = fieldnames(handles);

% Find all the elements whos font size we want to change:
textLocs = strfind(names, 'text_');
inputLocs = strfind(names, 'input_');
buttonLocs = strfind(names, 'button_');
panelLocs = strfind(names, 'panel_');
toggleLocs = strfind(names, 'toggle_');
popupLocs = strfind(names, 'popupmenu_');
iterLocs = strfind(names, 'iter_');
editLocs = strfind(names, 'edit_');

% Combine all the locations of elements whose font we wish to change
allLocs = strcat(textLocs, inputLocs, buttonLocs, panelLocs, toggleLocs, ...
    popupLocs, iterLocs, editLocs);

% Deal with variable maximum for input and noninput-type handles
for fieldCounter = 1:length(textLocs)
    if ( ~isempty(allLocs{fieldCounter}) )
        % Access the field values dynamically using the .( ) call
        currFS = get(handles.(names{fieldCounter}), 'FontSize');
        newFS = currFS + change;
        set(handles.(names{fieldCounter}), 'FontSize', newFS)
    end
end

% Change font sizes of the plots:
set(handles.fig_sol,  'fontsize', newFontsize)
set(handles.fig_norm, 'fontsize', newFontsize)

% Store the new font size:
handles.fontsizePanels = newFontsize;

end