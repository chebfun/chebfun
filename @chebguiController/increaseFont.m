function handles = increaseFont(handles)
% hObject    handle to button_fontinc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get a list of all the names in the handles
names = fieldnames(handles);
textLocs = strfind(names, 'text_');
inputLocs = strfind(names, 'input_');
buttonLocs = strfind(names, 'button_');
panelLocs = strfind(names, 'panel_');
toggleLocs = strfind(names, 'toggle_');

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

% Combine all the locations of elements whose font we wish to change
allLocs = strcat(textLocs, inputLocs, buttonLocs, panelLocs, toggleLocs);

for fieldCounter = 1:length(textLocs)
    if ( ~isempty(allLocs{fieldCounter}) )
        % Access the field values dynamically using the .( ) call. We
        % access the font size of each element separately as they can have
        % different relative font sizes.
        currFontSize = get(handles.(names{fieldCounter}), 'FontSize');
        if ( currFontSize >= 30 )
            % Do nothing (global max)
        elseif ( isempty(inputLocs{fieldCounter}) && (currFontSize >= 16) )
            % Size of non input-type handles is further limited
            if ( ~flag )
                % Record this, so a difference doesn't develop
                fontsizediff = min(fontsizediff + 1, 13); % 13 = 30-16-1
                flag = true; % but only do this once!
            end
        else
            newFontSize = currFontSize + 1;
            % Update the fontsize, again using the dynamic access
            set(handles.(names{fieldCounter}), 'FontSize', newFontSize)
        end
    end
end

handles.fontsizediff = fontsizediff;

end