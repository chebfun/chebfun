function chebtest(varargin)
%CHEBTEST   Run Chebfun test suite.
%
%   CHEBTEST executes all of the m-files found in the top level folders of the
%   directory <chebfunroot>/tests/. These m-files should return a scalar,
%   vector, or matrix of logical values. A test is deemed to pass if all the
%   returned values are logical true. There is no functional output from
%   CHEBTEST, but the data is piped to the command window with fprintf.
%
%   CHEBTEST('DIR1', 'DIR2', ...) will run only those tests given as inputs,
%   i.e., those in <chebfunroot>/tests/DIR1, <chebfunroot>/tests/DIR2, and so
%   on. If the given folder does not exist, a warning is given and that folder
%   skipped.
%
%   CHEBTEST will time each the tests, and if a test passes. If any of the
%   values returned by a test m-file are logical false (or zero) then the test
%   is deemed to have 'failed'. If a test results in an error being thrown it is
%   reported as 'crashed'. Therefore, an output from CHEBTEST might take the
%   form:
%     Test #001: chebtech1/test_alias.m ...          passed in 0.0094s
%     Test #002: chebtech1/test_bary.m ...           FAILED
%     Test #003: chebtech1/test_cell2mat.m ...       CRASHED
%     ...
%
%   If all test from a particular directory pass, the total execution time is
%   returned:
%     All chebtech1 tests passed in 3.1096s.
%   If not, the number of failed tests is reported:
%     1 failed test in chebtech2 directory.
%
%   Similarly, if all tests pass in all directories the total execution time is
%   returned at the end of the test:
%     All tests passed in 6.213s.
%   If not, then a list of the failed directories is reported:
%     Tests failed/crashed in directory: tests/chebtech2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]:  Preferences.

% Find directory in which Chebfun was installed:
installDir = fileparts(which('chebtest'));

% Set path to the tests/ subdirectory:
testsDir = [installDir '/tests'];

if ( nargin > 0 ) 
    % We are given a list of directories to test:
    testDirNames = varargin;
else
    % if not, run all the ones we can find in the tests/ folder.
    testDirNames = dir(testsDir);       % List test directories.
    testDirNames = {testDirNames.name}; % Keep only the directory names.
    testDirNames = testDirNames(3:end); % Strip '.' and '..' directories.
end

% If testDirNames is empty here, something is wrong.
if ( isempty(testDirNames) )
    error('CHEBFUN:chebtest:TestsNotFound', ...
        ['Could not locate test directories. ' ...
         'Please check that Chebfun has been installed correctly.']);
end

% Store the current directory: (We will return here when we're done.)
currDir = pwd();

% Switch to the tests/ subdirectory:
cd(testsDir);

% Initialise storage for which directories pass and the time they take.
numDirs = length(testDirNames);
passDir = zeros(numDirs, 1);
timeDir = zeros(numDirs, 1);

% Loop over the test directories and run the tests in each one.
for k = 1:numDirs
    testDir = testDirNames{k};
    if ( exist(testDir, 'dir') )
        fprintf(['Running tests in ' testDir ':\n']);
        [passDir(k), timeDir(k)] = runTestsInDirectory(testDir);
    else
        warning('CHEBFUN:chebtest:DirNotFound', ...
            'Test directory ''%s'' not found. Skipping.', testDir);
    end
    fprintf('\n');
end

if ( all(passDir) )
    % All tests have passed. Yay!
    fprintf('All tests passed in %4.4fs\n', sum(timeDir));
else
    % There's been a failure. List the directories which failed.
    if ( sum(~passDir) == 1 )
        % A single directory:
        fprintf('Tests failed/crashed in directory:\n   tests/%s.\n\n', ...
            testDirNames{~passDir});
    else
        % Multiple directories:
        failedDirsStr = '';
        for k = find(~passDir)'
            % Convert the list to a string:
            failedDirsStr = [failedDirsStr '\n   tests/' ...
                testDirNames{k}]; %#ok<AGROW>
        end
        fprintf(['Tests failed in directories:'  failedDirsStr '.\n\n']);
    end
end

% Return to current directory and return:
cd(currDir);

end

function [passFile, timeFile] = runTestsInDirectory(testDir)
%RUNTESTSINDIRECTORY   Run all the tests in the given directory.
%   RUNTESTSINDIRECTORY(TESTDIR) will change the current working directory to
%   TESTDIR, locate all the *.m files within it using dir, and execute each of
%   these in turn (in alphabetical order). If any of the tests crash (i.e.,
%   throw an error), it will be caught in a try-catch statement, and the error
%   will not be rethrown.
%
%   If all the tests in TESTDIR pass, then the total execution time for this
%   directory is also printed to screen.

% Store the current directory: (We will return here when we're done.)
currDir = pwd();

% Switch to the test directory and get list of all the test files:
cd(testDir);
testFiles = dir('*.m');
testFiles = {testFiles.name};

% For making the output string align nicely:
maxLength = max(cellfun(@length, testFiles));

% Allocate pass and timing variables:
numFiles = numel(testFiles);
passFile = zeros(numFiles, 1);
timeFile = zeros(numFiles, 1);

% Attempt to run all of the tests:
try % Note, we try-catch as we've CD'd and really don't want to end up elsewhere
    % Loop over the test files:
    for k = 1:numFiles
        % Next file to test: (.m extension is removed).
        testFile = testFiles{k}(1:end-2);
        printTestInfo(testDir, testFile, k, maxLength);
        [passFile(k), timeFile(k), resultStr] = runTest(testFile);
        fprintf([resultStr '\n']);
    end
catch ME
    % We failed. Return to the starting directory and rethrow the error:
    cd(currDir)
    rethrow(ME)
end

% Yay!
if ( all(passFile) )
    fprintf('All %s tests passed in %4.4fs.\n', testDir, sum(timeFile));
else
    fprintf('%d failed test in %s directory.\n', sum(~passFile), testDir);
end

% Restore the current working directory and return:
cd(currDir);

% Collapse to scalar for output:
passFile = all(passFile);
timeFile = sum(timeFile);

end


function printTestInfo(testDir, testFile, k, maxLength)
%PRINTTESTINFO  Pretty print test info for a given test file.
%
% PRINTTESTINFO(TESTDIR, TESTFILE, K) will print the following to screen:
%   Test #k: TESTDIR/TESTFILE.m ...
% Note that it will not linebreak after this.
%
% PRINTTESTINFO(TESTDIR, TESTFILE, K, MAXLENGTH) is similar, but will print
% extra whitespace after the ellipsis so that the following text for each
% TESTFILE in TESTDIR is aligned.
%
% Furthermore, PRINTTESTINFO tests to see if the Matlab desktop is running and
% will support HTML. If it is, the displayed the displayed TESTDIR/TESTFILE.m
% is hyperlinked to open TESTFILE.m in the Matlab editor.

if ( usejava('jvm') && usejava('desktop') )
    % Use HTML links if Java is enabled. 
    link = ['<a href = "matlab: edit ''', which(testFile), '''">' testDir '/', ...
        testFile '.m</a>'];                        
else
    % Otherwise, use plaintext.
    link = [testDir '/' testFile, '.m'];
end

ws = repmat(' ', 1, maxLength - length(testFile) - 1); % Whitespace.
fprintf('  Test #%3.3d: %s ...%s', k, link, ws);       % Print to screen.

end


function [pass, time, resultStr] = runTest(testFile)
%RUNTEST Runs the given test file.
%   [PASS, TIME, RESULTSTR] = RUNTEST(TESTFILE) executes the file TESTFILE.
%   TESTFILE should return a vector or matrix of logical values. 
%
%   If each of these are logical true (or 1) then the test passes and RUNTEST
%   returns PASS = TRUE, the time TIME that TESTFILE took to execute, and a
%   string RESULTSTR = 'passed in TIMEs'.
%
%   If any of the entries returned by TESTFILE are logcial false (or 0) then
%   PASS = FALSE and RESULTSTR = 'FAILED'.
%
%   If executing TESTFILE crashes, this is caught in a try-catch statement, and
%   RUNTTEST returns PASS = FALSE and RESULTSTR = 'CRASHED'.

% Attempt to run the test;
try
    tstart = tic();
    pass = feval(testFile);
    pass = all(pass(:));
    time = toc(tstart);
        
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

