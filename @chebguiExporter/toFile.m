function toFile(exporter, guifile, fileName, pathName)
%TOFILE     Export a CHEBGUI to an .m-file

% Concatenate the pathName and the fileName to get the full path:
fullFileName = [pathName, fileName];

try
    % Open a stream to write to a file:
    fid = fopen(fullFileName, 'wt');
    
    % Extract the necessary info for export to an .m file from the GUIFILE
    % object:
    expInfo = exporter.exportInfo(guifile);
    
    % Write the header information:
    writeHeader(exporter, fid, fileName)
    
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
    
    % Close the file writing stream:
    fclose(fid);
catch ME
    % Make sure to tidy up first
    fclose(fid);
    
    % Rethrow the error:
    rethrow(ME)
end

end