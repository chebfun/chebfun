function varargout = test(pref)
%TEST   Test the FUNCHEB2 class.
%   FUNCHEB2.TEST() runs the tests for the FUNCHEB2 class in the test/funcheb2/
%   folder. FUNCHEB2.test(PREF) overrides the default FUNCHEB2 preferences with
%   those in given in the preference structure PREF. P = FUNCHEB2.TEST returns a
%   boolean vector where the kth entry is true if the kth test passed.

%% 
% Set up.

% Initialise the random number generator:
rng('default')

% Store the current directory (we will return here when we're done).
currdir = pwd;

% Navigate to the funcheb2 test directory and choose the files to run:
funcheb2Dir = fileparts(which('funcheb2'));         % Locate funcheb2.
loc = find(funcheb2Dir == filesep, 1, 'last');      % Go down a directory.
testDir = fullfile(funcheb2Dir(1:loc), 'tests', 'funcheb2'); % Test directory.
cd(testDir)                                         % Move to it.
files = dir('*.m');                                 % Find all the m-files.
files = {files.name};                               % Ony want the name.

% Load some preferences:
if ( nargin == 0 )
    pref = funcheb.pref();
end    

% For making the output string align nicely:
maxLength = max(cellfun(@length, files));

% Pre-allocate the pass and timing vectors:
n = numel(files);
pass = zeros(1, n);
t = zeros(n, 1);

% If java is not enabled we won't display html links:
javaCheck = usejava('jvm') & usejava('desktop');

%% 
% Run the tests!

try % If we fail, we want to return to start directory.
    
    for k = 1:numel(files) % Loop over the test files:
        
        % Next file to test (.m extension is removed)
        file = files{k}(1:end-2);
        
        % Display information:
        if ( javaCheck )
            % Pretty:
            link = ['<a href = "matlab: edit ''', which(file), ...
                '''">funcheb2/', file '.m</a>'];                        
        else
            % Boring
            link = ['funcheb2/', file, '.m'];
        end
        ws = repmat(' ', 1, maxLength - length(file) - 1); % Whitespace.
        fprintf('  Test #%3.3d: %s ...%s', k, link, ws);   % Print to screen.
        
        try % If we fail, we want to report a crash.
            
            % Call the test:
            tic;
            pass(k) = all(feval(file, pref));
            t(k) = toc;
            
            % Did we pass?
            if ( pass(k) )
                % Hurrah!
                fprintf('passed in %2.4fs', t(k));
            else
                % :(
                fprintf('FAILED');
            end
            
        catch ME %#ok<NASGU>
            
            % We crashed. This is bad.
            fprintf('CRASHED');
            % But we _don't_ want to throw an error.
            %rethrow(ME)
            
        end
        
        % The next test goes on the next line:
        fprintf('\n');
        
    end % (for k = 1:numel(files))
    
    % Change back to the starting directory:
    cd(currdir)
    
catch ME
    
    % Change back to the starting directory:
    cd(currdir)
    
    % Rethrow the local error message we recieved:
    rethrow(ME)
    
end

% Yay!
if ( all(pass) )
    fprintf('All tests passed in %4.4fs\n', sum(t));
end

% Give an output argument if one was requested:
if ( nargout > 0 )
    varargout{1} = pass;
end

end
