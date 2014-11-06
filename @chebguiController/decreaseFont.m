function handles = decreaseFont(handles)
names = fieldnames(handles);
textLocs = strfind(names, 'text_');
inputLocs = strfind(names, 'input_');
buttonLocs = strfind(names, 'button_');
panelLocs = strfind(names, 'panel_');
toggleLocs = strfind(names, 'toggle_');

% Combine all the locations of elements whose font we wish to change
allLocs = strcat(textLocs, inputLocs, buttonLocs, panelLocs, toggleLocs);

% Deal with variable maximum for input and noninput-type handles
fontsizediff = [];
flag = 0;
tmp = strfind(names, 'fontsizediff');

if ( any([tmp{:}]) )
    fontsizediff = handles.fontsizediff;
end

if ( isempty(fontsizediff) )
    fontsizediff = 0;
end

for fieldCounter = 1:length(textLocs)
    if ( ~isempty(allLocs{fieldCounter}) )
        % Access the field values dynamically using the .( ) call
        currFontSize = get(handles.(names{fieldCounter}), 'FontSize');
        newFontSize = currFontSize;
        if ( isempty(inputLocs{fieldCounter}) )
            % Deal with non input-type handles
            if ( flag )
                % We only want update the diff once, so continue
                continue
            elseif ( fontsizediff > 0 )
                % Update the diff and set flag
                fontsizediff = fontsizediff - 1;
                flag = 1;
                continue
            elseif ( currFontSize > 5 )
                newFontSize = currFontSize - 1;
            end
        elseif ( currFontSize > 5 )
            % Input-type handles are easier
            newFontSize = currFontSize - 1;
        end
        % Update the fontsize, again using the dynamic access
        set(handles.(names{fieldCounter}), 'FontSize', newFontSize)
    end
end

handles.fontsizediff = fontsizediff;
end