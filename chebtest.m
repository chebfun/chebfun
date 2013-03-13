function chebtest(varargin)
%CHEBTEST   Run Chebfun test suite.
%    [TODO]:  Add documentation.

% [TODO]:  Preferences.

% Save current RNG state and re-initialize.
rngState = rng();
rng('default')

% Directory in which Chebfun was installed.
installDir = fileparts(which('chebtest'));

% Path to the tests/ subdirectory.
testsSubdir = [installDir '/tests'];

% If we were given a list of the directories for which to run the tests, use
% those.  Otherwise, run all of the ones we can find.
testDirNames = [];
if ( nargin == 0 ) testDirNames = dir(testsSubdir); % List test directories.
    testDirNames = {testDirNames.name}; % Keep only the directory names.
    testDirNames = testDirNames(3:end); % Strip '.' and '..' directories.
else
    testDirNames = varargin;
end

% If testDirNames is empty here, something is wrong.
if ( isempty(testDirNames) )
    error('CHEBFUN:chebtest:TestsNotFound', 'Could not locate test directories.  Please check that Chebfun has been installed correctly.');
end

% Store the current directory.  (We will return here when we're done.);
currDir = pwd();

% Switch to the tests/ subdirectory.
cd(testsSubdir);

% Loop over the test directories and run the tests in each one.
for ( n = 1:1:length(testDirNames) )
    testDir = testDirNames{n};
    if ( exist(testDir, 'dir') )
        fprintf(['Running tests in ' testDir ':\n']);
        runTestsInDirectory(testDir);
    else
        warning(['Test directory ' testDir ' not found.  Skipping.']);
    end
    fprintf('\n');
end

% Restore RNG state and current directory and return.
rng(rngState);
cd(currDir);

end


% [TODO]:  Add documentation.
function [pass, times] = runTestsInDirectory(testDir)

% Store the current directory.  (We will return here when we're done.);
currDir = pwd();

% Switch to the test directory and get list of all the test files.
cd(testDir);
testFiles = dir('*.m');
testFiles = {testFiles.name};

% For making the output string align nicely:
maxLength = max(cellfun(@length, testFiles));

% Allocate pass and timing variables.
n = numel(testFiles);
pass = zeros(1, n);
times = zeros(n, 1);

% Attempt to run all of the tests.
try
    % Loop over the test files.
    for ( k = 1:numel(testFiles) )
        % Next file to test (.m extension is removed).
        file = testFiles{k}(1:end-2);
        printTestInfo(testDir, file, k, maxLength);
        [pass(k), times(k), resultStr] = runTest(testDir, file);
        fprintf([resultStr '\n']);
    end
catch ME
    % We failed.  Return to the starting directory and rethrow the error.
    cd(currDir)
    rethrow(ME)
end

% Yay!
if ( all(pass) )
    fprintf('All tests passed in %4.4fs\n', sum(times));
end

% Restore the current working directory and return.
cd(currDir);

end


% [TODO]:  Add documentation.
function printTestInfo(testDir, file, k, maxLength)

% Use HTML links if Java is enabled.  Otherwise, use plaintext.
if ( usejava('jvm') && usejava('desktop') )
    link = ['<a href = "matlab: edit ''', which(file), ...
            '''">' testDir '/', file '.m</a>'];                        
else
    link = [testDir '/' file, '.m'];
end

ws = repmat(' ', 1, maxLength - length(file) - 1); % Whitespace.
fprintf('  Test #%3.3d: %s ...%s', k, link, ws);   % Print to screen.

end


% [TODO]:  Add documentation.
function [pass, time, resultStr] = runTest(testDir, file)

% Attempt to run the test.
try
    tic();
    pass = feval(file);
    pass = all(pass(:));
    time = toc();
        
    % Did we pass?
    if ( pass )
        % Hurrah!
        resultStr = sprintf('passed in %2.4fs', time);
    else
        % :(
        resultStr = 'FAILED';
    end

catch ME %#ok<NASGU>
    % We crashed. This is bad.
    pass = false;
    time = NaN;
    resultStr = 'CRASHED';

    % But we _don't_ want to throw an error.
    %rethrow(ME)
end

end

