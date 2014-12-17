function pass = test_toFilePDE(pref)
%TEST_TOFILEPDE   Test exporting all PDE demos to an m-file.
%
% This test only checks whether nothing breaks, it does not try to solve the
% problems.

% This test won't work if we can't write to the current directory.
[ignored, attr] = fileattrib(pwd);
if ( ~attr.UserWrite )
    warning('CHEBFUN:tests:chebgui:test_toFilePDE:perms', ...
        'Cannot write to chebfunroot/tests/chebgui/.  Bypassing test.');
    pass = true;
    return;
end

% Find the folders which demos are stored in. The chebguiDemos folder lives in
% the trunk folder, find the path of the Chebfun trunk.
trunkPath = chebfunroot();

% Append directory information
bvppath = fullfile(trunkPath, 'chebguiDemos', 'pdedemos');

% Set up file to exporting to
tempPath = fullfile(trunkPath, 'tests', 'chebgui');
tempFileName = 'tempExporterTest.m';

% Obtain directory information:
D = dir(bvppath);

% Create a BVP exporter object:
exporter = chebguiExporter.constructor('pde');

% Wrap the test in a try-catch statement to ensure we don't files hanging around
% in the test folder:
try
    
    % Loop through the demos and export
    for demoCounter = 3:length(D)
        
        % The path to the .guifile
        demoPath = fullfile(bvppath, D(demoCounter).name);
        
        % Load demo
        cg = chebgui.demo2chebgui(demoPath);
        
        % Export the demo!
        toFile(exporter, cg, tempFileName, tempPath)
    end
    
    % Got here without crashing == success!
    pass = 1;

catch
    % Something went wrong == fail :(
    pass = 0;
end

% Delete the temporary file we wrote to:
delete(fullfile(tempPath, tempFileName))


end
