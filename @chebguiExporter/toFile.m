function toFile(exporter, guifile)

[filename, pathname] = uiputfile( ...
    {'*.m', 'M-files (*.m)'; '*.*',  'All Files (*.*)'}, ...
    'Save as', exporter.defaultFileName);

fullFileName = [pathname, filename];

fid = fopen(fullFileName, 'wt');

% Extract the necessary info for export to an .m file from the GUIFILE object:
expInfo = exportInfo(exporter, guifile);

writeHeader(exporter, fid, filename)

if ( filename ~= 0 )     % User did not press cancel
    % Print description of the problem:
    exporter.printDescription(fid, expInfo)
    
    % Print the problem set-up:
    exporter.printSetup(fid, expInfo, guifile)
    
    % Print the options set-up:
    exporter.printOptions(fid, expInfo)
    
    % Print the solution step:
    exporter.printSolver(fid, expInfo)
    
    % Print the post-solution process:
    exporter.printPostSolver(fid, expInfo)
    
    % Open the new file in the editor
    open(fullFileName)
end


fclose(fid);


end