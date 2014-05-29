function toFile(e, guifile)

[filename, pathname] = uiputfile( ...
    {'*.m', 'M-files (*.m)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Save as', e.defaultFileName);

if ( filename ~= 0 )     % User did not press cancel
    try
        chebgui2mfile(e, guifile, pathname, filename)
        
    catch ME
        rethrow(ME)
        error('Chebgui:Export', ...
            ['Error in exporting to .m file. Please make ' ...
            'sure there are no syntax errors.']);
    end
    
    % Open the new file in the editor
    open([pathname, filename])
end



end