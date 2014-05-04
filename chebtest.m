function varargout = chebtest(varargin)
%CHEBTEST   Run Chebfun test suite.
%
%   RESULTS = CHEBTEST executes all of the m-files found in the top level
%   folders of the directory <chebfunroot>/tests/. These m-files should return
%   a scalar, vector, or matrix of logical values. A test is deemed to pass if
%   all the returned values are logical true. If the output argument RESULTS is
%   an N-by-3 cell-array with one row per test and columns:
%     (1) directory name, a string.     Example: 'chebtech2'
%     (2) test name, a string.          Example: 'test_chebpts'
%     (3) duration OR error flag. A value >0 indicates success and is equal to
%         the execution time of the test in seconds, while -1 == failure and
%         -2 == crash.
%
%   CHEBTEST('DIR1', 'DIR2', ...) will run only those tests given as inputs,
%   i.e., those in <chebfunroot>/tests/DIR1, <chebfunroot>/tests/DIR2, and so
%   on. If the given folder does not exist, a warning is given and that folder
%   skipped.
%
%   CHEBTEST(..., '--log') will write the results of the tests to the file
%   "chebtest-YYYYMMDDHHMMSS.log" in the current directory, with
%   "YYYYMMDDHHMMSS" a numeric datetime.
%
%   CHEBTEST will time each of the tests, and print the time for each test
%   that passes. If any of the values returned by a test m-file are logical
%   false (or zero) then the test is deemed to have 'failed'. If a test
%   results in an error being thrown it is reported as 'crashed'. Therefore,
%   an output from CHEBTEST might take the form:
%   | Test #001: chebtech1/test_alias.m ...          passed in 0.0094s
%   | Test #002: chebtech1/test_bary.m ...           FAILED
%   | Test #003: chebtech1/test_cell2mat.m ...       CRASHED
%   | ...
%
%   If all test from a particular directory pass, the total execution time is
%   returned:
%   | All chebtech1 tests passed in 3.1096s.
%   If not, the number of failed tests is reported:
%   | 1 failed test in chebtech2 directory.
%
%   Similarly, if all tests pass in all directories the total execution time is
%   returned at the end of the test:
%   | All tests passed in 6.213s.
%   If not, then a list of the failed directories is reported:
%   | The following tests failed or crashed:
%   |    chebtech2/test_chebpts
%   |    adchebfun/test_erf

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Preferences.

% Find directory in which Chebfun was installed:
installDir = fileparts(which('chebtest'));

% Set path to the tests/ subdirectory:
testsDir = [installDir '/tests'];

% Set the flag for logging results.
writeLog = false;
if ( nargin > 0 )
    % We are given a list of directories to test:
    if ( strcmp(varargin{end}, '--log') )
        writeLog = true;
        varargin = {varargin{1:end-1}};
    end
end

% Sort out other input arguments, the directories to test.
if ( length(varargin) )
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
allResults = cell(0,3);

% Loop over the test directories and run the tests in each one.
for k = 1:numDirs
    testDir = testDirNames{k};
    if ( exist(testDir, 'dir') )
        fprintf(['Running tests in ' testDir ':\n']);
        nextResults = runTestsInDirectory(testDir);
        allResults = [allResults; nextResults];
    else
        warning('CHEBFUN:chebtest:DirNotFound', ...
            'Test directory ''%s'' not found. Skipping.', testDir);
    end
    fprintf('\n');
end

durations = [allResults{:,3}];
if ( all(durations > 0) )
    % All tests have passed. Yay!
    fprintf('All tests passed in %4.4fs\n', sum(durations));

else
    % There's been a failure. 

    indx = find(durations < 0);
    fprintf('The following tests failed or crashed:\n');
    for k = indx(:)'
        fprintf('   %s/%s\n', allResults{k,1:2})
    end

end

% Return to current directory.
cd(currDir);

% Write the log if requested to.
if ( writeLog )
    filename = ['chebtest-' datestr(now, 'yyyymmddHHMMSS') '.log'];
    writeToLog(filename, allResults);
end

% Return.
if ( nargout > 0 )
    varargout{1} = allResults;
end

end

function testResults = runTestsInDirectory(testDir)
%RUNTESTSINDIRECTORY   Run all the tests in the given directory.
%   RUNTESTSINDIRECTORY(TESTDIR) will change the current working directory to
%   TESTDIR, locate all the *.m files within it using DIR, and execute each of
%   these in turn (in alphabetical order). If any of the tests crash (i.e.,
%   throw an error), it will be caught in a try-catch statement, and the error
%   will not be rethrown.
%
%   If all the tests in TESTDIR pass, then the total execution time for this
%   directory is also printed to screen.
%
%   RUNTESTSINDIRECTORY returns a cell-array with values the same as described
%   in the documentation for WRITETOLOG below.

% Store the current directory: (We will return here when we're done.)
currDir = pwd();

% Switch to the test directory and get list of all the test files:
cd(testDir);
testFiles = dir('*.m');
testFiles = {testFiles.name};

% For making the output string align nicely:
maxLength = max(cellfun(@length, testFiles));

% Allocate pass and timing variables:
numFiles  = numel(testFiles);
durations = zeros(numFiles, 1);
errorMessages = {'FAILED', 'CRASHED'};

% Attempt to run all of the tests:
try % Note, we try-catch as we've CD'd and really don't want to end up elsewhere
    % Loop over the test files:
    for k = 1:numFiles
        % Next file to test: (.m extension is removed).
        testFile = testFiles{k}(1:end-2);
        printTestInfo(testDir, testFile, k, maxLength);
        durations(k) = runTest(testFile);
        if durations(k) > 0,
            % Success message.
            message = sprintf('passed in %.4fs', durations(k));
        else
            % Error flags are negative indices.
            message = errorMessages{-durations(k)};
        end
        fprintf([message '\n']);
    end
catch ME
    % We failed. Return to the starting directory and rethrow the error:
    cd(currDir)
    rethrow(ME)
end

% Yay!
if ( all(durations > 0) )
    fprintf('All %s tests passed in %4.4fs.\n', testDir, sum(durations));
else
    fprintf('%d failed test in %s directory.\n', sum(durations < 0), testDir);
end

% Restore the current working directory and return:
cd(currDir);

% Dump all the test data into a cell array to pass back to the CHEBTEST
% function. This data is what is written to a .CSV file if logging is on.
testResults = cell(numFiles,3);
for k = 1:numFiles
    testResults{k,1} = testDir;               % directory name
    testResults{k,2} = testFiles{k}(1:end-2); % test file name
    testResults{k,3} = durations(k);          % duration / error flag
end

end


function printTestInfo(testDir, testFile, k, maxLength)
%PRINTTESTINFO  Pretty print test info for a given test file.
%
%   PRINTTESTINFO(TESTDIR, TESTFILE, K) will print the following to screen:
%     Test #k: TESTDIR/TESTFILE.m ...
%   Note that it will not linebreak after this.
%
%   PRINTTESTINFO(TESTDIR, TESTFILE, K, MAXLENGTH) is similar, but will print
%   extra whitespace after the ellipsis so that the following text for each
%   TESTFILE in TESTDIR is aligned.
%
%   Furthermore, PRINTTESTINFO tests to see if the Matlab desktop is running and
%   will support HTML. If it is, the displayed the displayed TESTDIR/TESTFILE.m
%   is hyperlinked to open TESTFILE.m in the Matlab editor.

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


function duration = runTest(testFile)
%RUNTEST Runs the given test file.
%   DURATION = RUNTEST(TESTFILE) executes the file TESTFILE.
%   TESTFILE should return a vector or matrix of logical values. 
%
%   If each of these are logical true (or 1) then the test passes and RUNTEST
%   returns DURATION as a double > 0 corresponding to the time it took the test
%   to execute in seconds.
%
%   If any of the entries returned by TESTFILE are logcial false (or 0), then
%   DURATION = -1.
%
%   If executing TESTFILE crashes, this is caught in a try-catch statement, and
%   RUNTTEST returns DURATION = -2.

% Attempt to run the test:
try
    tstart = tic();
    pass = feval(testFile);
    duration = toc(tstart);

    pass = all(pass(:));
    if ~pass,
        % The test failed, so return FAILED flag.
        duration = -1;
    end

catch ME %#ok<NASGU>
    % We crashed. This is bad. Return CRASHED flag.
    duration = -2;

    % But we _don't_ want to throw an error.
    %rethrow(ME)
end

end


function writeToLog(filename, data)
%WRITETOLOG   Writes the contents of a cell-array to a file as CSV.
%
%   WRITETOLOG(FN,DATA) writes the contents of the N-by-3 cell-array DATA to
%   the file FN. The cell-array DATA must have three columns:
%     (1) Directory name, a string.     Example: 'chebtech2'
%     (2) Test name, a string.          Example: 'test_chebpts'
%     (3) Duration in seconds OR error flag. A value >0 indicates success,
%         while -1 == failure and -2 == crash.

columnTitles = {'dir_name', 'test_name', 'duration'};
data = data';

fid = fopen(filename, 'w');
fprintf(fid, '%s,%s,%s\n', columnTitles{:});
fprintf(fid, '%s,%s,%f\n', data{:});
fclose(fid);

end
