function printOptions(fid, expInfo)
%PRINTOPTIONS   Print problem options when they are exported.
%
% Calling sequence:
%   PRINTOPTIONS(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract info from the expInfo struct:
K = expInfo.K;
discretization = expInfo.discretization;

% Set up preferences
fprintf(fid, '\n%%%% Setup preferences for solving the problem.\n');
fprintf(fid, '%% Create a CHEBOPPREF object for passing preferences.\n');
fprintf(fid, '%% (See ''help cheboppref'' for more possible options.)\n');
fprintf(fid, 'options = cheboppref();\n');

% Option for discretization:
fprintf(fid, '\n%% Specify the discretization to use. Possible options are:\n');
fprintf(fid, '%%  ''values'' (default)\n');
fprintf(fid, '%%  ''coeffs''\n');
fprintf(fid, '%%  A function handle (see ''help cheboppref'' for details).\n');
fprintf(fid, 'options.discretization = ''%s'';\n', discretization);

% Specify number of eigenvalues to compute:
fprintf(fid, '\n%% Number of eigenvalue and eigenmodes to compute.\n');
fprintf(fid, 'k = %s;\n', K);

end
