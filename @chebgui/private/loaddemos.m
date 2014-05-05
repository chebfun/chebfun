function cg = loaddemos(guifile,guifilepath) %#ok<INUSL>
% Load a demo from a .guifile to a chebgui object

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Import from the given file and evaluate to fill the workspace
fid = fopen(guifilepath);

% Make sure the file exists
if fid == -1
  error('CHEBGUI:noload','Unable to open demo file: %s.',guifilepath)
end

% Load the data to the workspace
inputEnded = 0;
while ~inputEnded
    tline = fgetl(fid);
    if isempty(strfind(tline,'=')), continue, end % Don't eval names and demotypes
    if ~ischar(tline), break, end
    eval(tline);
    inputEnded = feof(fid);
end
fclose(fid);

% timedomain is entered as 't' in .guifiles. Sort this out:
if exist('t','var')
    timedomain = t; %#ok<NASGU>
    clear t
end

% Clear these variables that we're finished with
clear inputEnded fid tline ans

% Load all the current workspace vars into the chebgui object using SET
vars = who;
cg = chebgui('type','bvp');
for k = 1:numel(vars)
    if any(strcmp(vars{k},{'guifile','guifilepath'})), continue, end
    try
        cg = set(cg,vars{k},eval(vars{k}));
    catch ME %#ok<NASGU>
        warning('CHEBGUI:loaddemos:unknown',...
            [vars{k} ' is an unknown CHEBGUI property. Ignoring.']);
    end
end

