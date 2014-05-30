function chebgui2mfile(exporter, guifile, fid, expInfo)
%CHEBGUI2MFILE      Export a BVP an .m-file

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO:  Documentation.

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

end