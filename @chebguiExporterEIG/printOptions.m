function printOptions(fid, expInfo)
%PRINTOPTIONS   Print problem options when they are exported.
%
% Calling sequence:
%   PRINTOPTIONS(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract info from the expInfo struct:
K = expInfo.K;
discretization = expInfo.discretization;

% Set up preferences
fprintf(fid, '\n%%%% Setup preferences for solving the problem.');
fprintf(fid, '\n%% Create a CHEBOPPREF object for passing preferences: \n');
fprintf(fid, 'options = cheboppref();\n');

% Option for discretization:
fprintf(fid, ['\n%% Option for discretization (either ''collocation'' ' ...
    'or ''ultraspherical'').\n']);
if ( isa(discretization(), 'colloc') )
    discString = 'collocation';
else
    discString = 'ultraspherical';
end
fprintf(fid, 'options.discretization = ''%s'';\n', discString);
fprintf(fid, '\n%% Number of eigenvalue and eigenmodes to compute.\n');
fprintf(fid, 'k = %s;\n', K);

end