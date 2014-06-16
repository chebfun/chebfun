function pass = test_toFileEIG(pref)
%TEST_TOFILEEIG     Test exporting all EIG demos to an .m-file.
%
% This test only checks whether nothing breaks, it does not try to solve the
% problems.

% Find the folders which demos are stored in. The chebguiDemos folder lives in
% the trunk folder, find the path of the Chebfun trunk.
trunkPath = chebfunroot();

% Append directory information
bvppath = fullfile(trunkPath, 'chebguiDemos', 'eigdemos');

% Set up file to exporting to
tempPath = fullfile(trunkPath, 'tests', 'chebgui');
tempFileName = 'tempExporterTest.m';

% Obtain directory information:
D = dir(bvppath);

% Create a BVP exporter object:
exporter = chebguiExporter.constructor('eig');

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