function handles = changeFontsize(handles, change)

currentChange = handles.fontsizeChanges;
newChange = currentChange + change;

% Allow max 5 increases or 5 decreases
maxChange = 5;
minChange = -5;

if ( (newChange > maxChange) || (newChange < minChange) )
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

% Combine all the locations of elements whose font we wish to change
allLocs = strcat(textLocs, inputLocs, buttonLocs, panelLocs, toggleLocs);

% Deal with variable maximum for input and noninput-type handles
for fieldCounter = 1:length(textLocs)
    if ( ~isempty(allLocs{fieldCounter}) )
        % Access the field values dynamically using the .( ) call
        currFontSize = get(handles.(names{fieldCounter}), 'FontSize');
        newFontSize = currFontSize + change;
        set(handles.(names{fieldCounter}), 'FontSize', newFontSize)
    end
end

handles.fontsizeChanges = newChange;

end