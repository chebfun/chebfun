function cg = demo2chebgui(demoPath)
%DEMO2CHEBGUI   Load a demo stored in a .guifile to a CHEBGUI object
%   CG = CHEBGUI.DEMO2CHEBGUI(DEMOPATH) convert the .guifile stored on DEMOPATH
%   to a CHEBGUI object.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Import from the given file and evaluate to fill the workspace
fid = fopen(demoPath);

% Make sure the file exists:
if ( fid == -1 )
    error('CHEBFUN:CHEBGUI:demo2chebgui:noload', ...
        'Unable to open demo file: %s.', demoPath)
end

% Load the data to the workspace
inputEnded = 0;
while ( ~inputEnded )
    tline = fgetl(fid);

    % Don't eval names and demotypes. Also, don't eval lines that start at #, as
    % those are problem description lines needed for testing.
    if ( isempty(tline) || isempty(strfind(tline, '=')) || tline(1) == '#' )
        continue
    end

    if ( ~ischar(tline) )
        break
    end

    eval(tline);
    inputEnded = feof(fid);
end

% Close the file stream:
fclose(fid);

% timedomain is entered as 't' in .guifiles. Sort this out:
if ( exist('t', 'var') )
    timedomain = t; %#ok<NASGU>
    clear t
end

% Clear these variables, that we're finished with
clear inputEnded fid tline ans

% Load all the current workspace vars into the chebgui object using SET:
vars = who;
cg = chebgui('type', 'bvp');
for k = 1:numel(vars)
    if ( strcmp(vars{k}, 'demoPath') )
        continue
    end
    try
        cg = set(cg, vars{k}, eval(vars{k}));
    catch ME %#ok<NASGU>
        warning('CHEBFUN:CHEBGUI:loaddemos:unknown',...
            [vars{k} ' is an unknown CHEBGUI property. Ignoring.']);
    end
end

end
