function pass = test_toFileEIG(pref)
%TEST_TOFILE    Test exporting all BVP demos to an .m-file.
%
% This test only checks whether nothing breaks, it does not try to solve the
% problems.
% Find the folders which demos are stored in. The chebguiDemos folder lives in
% the trunk folder, find the path of the Chebfun trunk.
trunkPath = fileparts(which('chebguiWindow'));

% Append directory information
bvppath = fullfile(trunkPath, 'chebguiDemos', 'eigdemos');

% Set up file to exporting to
tempPath = fullfile(trunkPath, 'tests', 'chebguiExporter');
tempFileName = 'tempExporterTest.m';

% Obtain directory information:
D = dir(bvppath);

% Create a BVP exporter object:
exporter = chebguiExporter.constructor('eig');

% Loop through the demos and export
for demoCounter = 3:length(D)
    demoPath = fullfile(bvppath, D(demoCounter).name);
    
    % Load demo
    cg = chebgui.demo2chebgui(demoPath);
        
    % Export the demo!
    toFile(exporter, cg, tempFileName, tempPath)
end

% Delete the temporary file we wrote to:
delete(fullfile(tempPath, tempFileName))

% Got here without crashing == success!
pass = 1;

end