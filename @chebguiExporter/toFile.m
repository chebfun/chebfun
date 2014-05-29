function toFile(e, guifile)

[filename, pathname] = uiputfile( ...
    {'*.m', 'M-files (*.m)'; '*.*',  'All Files (*.*)'}, ...
    'Save as', e.defaultFileName);

fullFileName = [pathname, filename];

fid = fopen(fullFileName, 'wt');

exportInfo(e, fid, filename)

if ( filename ~= 0 )     % User did not press cancel
    chebgui2mfile(e, guifile, fid)
            
    % Open the new file in the editor
    open(fullFileName)
end


fclose(fid);


end