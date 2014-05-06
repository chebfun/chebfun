function updateChebfun(varargin)
%UPDATECHEBFUN   Update Chebun source files.
%   UPDATECHEBFUN() updates the Chebfun source files with the latest stable
%   relase from the Chebfun Github repository https://github.com/chebfun/chebfun
%   WARNING: ALL FOLDERS AND FILES IN THE $CHEBFUNROOT DIRECTORY WILL BE REMOVED
%   IN THE PROCESS.
%
%   UPDATECHEBFUN('devel') updates to the latest development build. Rather than
%   use this, it is recommended that you clone a copy copy of the repository.
%
%   UPDATECHEBFUN(..., '--backup') makes a backup .zip file of the $chebfunroot
%   folder before removing and replacnig it.
%
%   UPDATECHEBFUN(..., '--force') skips the warning dialogues before updating.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: This won't work until the Git repo becomes public.

version = 'release';  % 'release' or 'development'?
force = false;        % Don't ask for confirmation.
backup = false;       % Make a .zip of the current $chebfunroot directory.

% Parse inputs:
for k = 1:numel(varargin)
    if ( any(strcmpi(varargin{k}, {'release', 'stable'})) )
        version = 'stable';
    elseif ( any(strcmpi(varargin{k}, {'devel', 'development', 'nightly'})) )
        version = 'development';
    elseif ( any(strcmpi(varargin{k}, {'--force', '--f'})) )
        force = true;
    else
        error('CHEBFUN:updateChebfun:unknown', ...
            'Unknown input %s.', varargin{k});
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check they _really_ want to do this!
if ( ~force )
    
    fprintf(2, 'Warning: Updating will permenently remove your entire Chebfun root directory and all contents.\n');
    yesno = '';
    while ( ~any(strcmp(yesno, {'yes', 'no'}))  )
        fprintf(2, 'Proceed (yes / no)? ');
        yesno = input('', 's');
        if ( strcmpi(yesno, 'no') ), return, end
    end
    
    % Check again!
    yesno = '';
    while ( ~any(strcmp(yesno, {'yes', 'no'}))  )
        fprintf(2, 'Are you sure (yes / no)? ');
        yesno = input('', 's');
        if ( strcmpi(yesno, 'no') ), return, end
    end
    
    % Check for a Git repo:
    if ( exist('.git', 'dir') )
        fprintf(2, 'Warning: Your Chebfun root directory appears to be part of a Git repo.\n');
        yesno = '';
        while ( ~any(strcmp(yesno, {'yes', 'no'}))  )
            fprintf(2, 'Do you still want to continue (yes / no)? ');
            yesno = input('', 's');
            if ( strcmpi(yesno, 'no') ), return, end
        end
    end
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try % Wrap everything in a try-catch statement to avoid a disaster.
    
    % Remember the current directory:
    currDir = pwd;
    % Find directory in which Chebfun was installed:
    installDir = fileparts(which('chebtest'));
    % Navigate there:
    cd(installDir)
    
    if ( backup ) %#ok<*UNRCH>
        % Make a back up, so that we can recover things if needs be.
        filename = ['chebfun_backup_', datestr(now, 'yyyymmddHHMMSS'), '.zip']; 
        warnState = warning('off', 'MATLAB:zip:archiveName');
        zip(filename, '*', pwd);
        warning(warnState);
    end

    % Navigate down a directory:
    cd ..
    % Remove the Chebfun directory:
    disp(['Removing ', installDir, '...']);
%     rmdir(installDir, 's'); % TODO: Replace this line.
    
    % Get the requested version of Chebfun:
    githuburl = 'https://github.com/chebfun/chebfun/';
    if ( strcmp(version, 'release') )
        disp('Downloading .zip file of latest stable release... ')
        pause(.5)
%         unzip([githuburl, 'archive/master.zip'], installDir);
    else
        disp('Downloading .zip file of latest development version... ')
        pause(.5)
%         unzip([githuburl, 'archive/development.zip'], installDir);
    end
    disp('Extracting .zip file...')
    
    % Save the path:
    addpath(installDir)
    savepath
    % Return to the starting directory:
    cd(currDir)

    % Tidy up:
    close all
    clear all
    clear classes
    
    % Display outgoing info:
    disp('Installation complete!')
    
    % Suggest they run CHEBTEST():
    if ( usejava('jvm') && usejava('desktop') )
        str = '<a href="matlab: chebtest">chebtest</a>';
    else
        str = 'chebtest';
    end
    disp(['We recommend you run ' str ' to ensure everything is working.'])
    
catch ME
    
    % Something went wrong. Throw an error:
    disp('Update failed.');
    rethrow(ME)
    
    % Not that even if the $chebfunroot folder was removed, it is recovered by
    % the try-catch statement.
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


    

