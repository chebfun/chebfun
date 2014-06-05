function pass = test_toFile(pref)
%TEST_TOFILE    Test exporting a selection of problems to .m files
%
% This test only checks whether nothing breaks, it does not try to solve the
% problems.
% Find the folders which demos are stored in. The chebguiDemos folder lives in
% the trunk folder, find the path of the Chebfun trunk.
trunkPath = fileparts(which('chebguiWindow'));

% Append directory information
bvppath = fullfile(trunkPath, 'chebguiDemos', 'bvpdemos');
eigpath = fullfile(trunkPath, 'chebguiDemos', 'eigdemos');
pdepath = fullfile(trunkPath, 'chebguiDemos', 'pdedemos');


% Set up file to exporting to
tempPath = fullfile(trunkPath, 'tests', 'chebgui');
tempFileName = 'tempExporterTest.m';

demos = {'interior_layers_nonlinear.guifile'; 'system_coupled_bcs.guifile'; ...
    'mathieu.guifile'; 'harmonic_system.guifile'; ...
    'allen_cahn.guifile'; 'diffusion_three_chem.guifile'};
types = {'bvp'; 'bvp'; 'eig'; 'eig'; 'pde'; 'pde'};

try
    % We wrap this whole test in a try-catch so that we don't leave files
    % hanging around in the test folder.

    % Loop through the problems we want to export:
    for demoCounter = 1:numel(demos)

        % Create an exporter object of the correct type:
        exporter = chebguiExporter.constructor(types{demoCounter});

        % Full path to demo
        if demoCounter <=2 
            demoPath = fullfile(bvppath, demos{demoCounter});
        elseif demoCounter <= 4
            demoPath = fullfile(eigpath, demos{demoCounter});
        else
            demoPath = fullfile(pdepath, demos{demoCounter});
        end
        % Load demo
        cg = chebgui.demo2chebgui(demoPath);

        % Export the demo!
        toFile(exporter, cg, tempFileName, tempPath)

    end
    
    % Got here without crashing == success!
    pass = 1;
    
catch
    
    % Something went wrong == fail.
    pass = false;
    
end

% Delete the temporary file we wrote to:
delete(fullfile(tempPath, tempFileName))

end