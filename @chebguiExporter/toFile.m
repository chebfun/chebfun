function toFile(exporter, guifile, fileName, pathName)
%TOFILE     Export a CHEBGUI to an .m-file
%
% Calling sequence:
%
%   TOFILE(EXPORTER, GUIFILE, FILENAME, PATHNAME)
%
% where
%
%   EXPORTER:   A CHEBGUIEXPORTER object (of BVP, EIG or PDE kind).
%   GUIFILE:    A CHEBGUI object.
%   FILENAME:   The name of the file we want to write to.
%   PATHNAME:   The path to the file we want to write to.
%
% This method will export the problem to the file 'PATHNAME/FILENAME'.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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