function chebgui2mfile(exporter, guifile, fid, expInfo)
%CHEBGUI2MFILE      Export a BVP an .m-file

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO:  Documentation.

% Extract information from the EXPINFO struct
allVarNames = expInfo.allVarNames;
indVarNameSpace = expInfo.indVarNameSpace;

% Print description of the problem:
exporter.printDescription(fid, expInfo)

% Print the problem set-up:
exporter.printSetup(fid, expInfo, guifile)

% Set up preferences
fprintf(fid, '\n%% Setup preferences for solving the problem.\n');
fprintf(fid, 'options = cheboppref();\n');

% Option for tolerance
tolInput = guifile.tol;

if ( ~isempty(tolInput) )
    fprintf(fid, '\n%% Option for tolerance.\n');
    fprintf(fid, 'options.errTol = %s;\n', tolInput);
end

fprintf(fid, 'options.display = ''iter'';\n');

% Option for damping
dampedOnInput = guifile.options.damping;

fprintf(fid, '\n%% Option for damping.\n');
if ( strcmp(dampedOnInput, '1') )
    fprintf(fid, 'options.damped = ''on'';\n');
else
    fprintf(fid, 'options.damped = ''off'';\n');
end

fprintf(fid, '\n%% Option for discretization (either @colloc2 or @ultraS).\n');
fprintf(fid, 'options.discretization = @%s;\n', ...
    func2str(guifile.options.discretization));

% Option for plotting
plottingOnInput = guifile.options.plotting;

if ( ~strcmp(plottingOnInput, 'off') )
    fprintf(fid, ['\n%% Option for determining how long each Newton step ' ...
        'is shown.\n']);
    fprintf(fid, 'options.plotting = %s;\n', plottingOnInput);
else
    fprintf(fid, ['\n%% Option for determining how long each Newton step ' ...
        'is shown.\n']);
    fprintf(fid,'options.plotting = ''off'';\n');
end

fprintf(fid, ['\n%% Solve the problem using solvebvp (a routine which ' ...
    'offers the same\n%% functionality as nonlinear backslash, but with '...
    'more customizability).\n']);
fprintf(fid, 'u = solvebvp(N, rhs, options);\n');

fprintf(fid, '\n%% Create a plot of the solution.\n');

fprintf(fid, ['figure\nplot(u,''LineWidth'',2)\ntitle(''Final solution'')', ...
    ', xlabel(''%s'')'], indVarNameSpace);
if ( numel(allVarNames) == 1 )
    fprintf(fid, ', ylabel(''%s'')', allVarNames{:});
else
    leg = '';
    for k = 1:numel(allVarNames)-1
        leg = [leg '''' allVarNames{k} '''' ','];
    end
    leg = [leg '''' allVarNames{k+1} ''''];
    fprintf(fid, ', legend(%s)\n', leg);
end


end