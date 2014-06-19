function update(varargin)
%CHEBFUN.UPDATE   Update Chebun source files.
%   UPDATECHEBFUN() updates the Chebfun source files with the latest stable
%   release from the Github repository https://github.com/chebfun/chebfun/
%   WARNING: ALL FOLDERS AND FILES IN THE $CHEBFUNROOT/ DIRECTORY WILL BE
%   REMOVED IN THE PROCESS.
%
%   CHEBFUN.UPDATE('devel') updates to the latest development build. Rather than
%   use this, it is recommended that you clone a copy of the repository. See
%   http://chebfun.github.io/develop/ for details.
%
%   CHEBFUN.UPDATE(..., '--backup') makes a backup .zip file of the
%   $chebfunroot/ folder before removing and replacing it.
%
%   CHEBFUN.UPDATE(..., '--force') skips the warning dialogues before updating.
%
%   Note that this file will also remove all current variables from the
%   workspace, as well as compiled MATLAB and MEX-functions. (See help clear).
%
% See also CHEBFUNROOT().

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

version = 'release';  % 'release' or 'development'?
force = false;        % Don't ask for confirmation.
backup = false;       % Make a .zip of the current $chebfunroot/ directory.
gitHubURL = 'https://github.com/chebfun/chebfun/'; % GitHub URL

% Parse inputs:
for k = 1:numel(varargin)
    if ( any(strcmpi(varargin{k}, {'release', 'stable'})) )
        version = 'stable';
    elseif ( any(strcmpi(varargin{k}, ...
            {'devel', 'development', 'nightly', 'dev'})) )
        version = 'development';
    elseif ( any(strcmpi(varargin{k}, {'--force', '--f'})) )
        force = true;
    elseif ( any(strcmpi(varargin{k}, {'--backup'})) )
        backup = true;        
    else
        error('CHEBFUN:updateChebfun:unknown', ...
            'Unknown input %s.', varargin{k});
    end
end

% Remember the current directory:
startDir = pwd;
% Find directory in which Chebfun was installed:
installDir = chebfunroot();
% Navigate there:
cd(installDir)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check they _really_ want to do this!
if ( ~force )
    
    fprintf(2, ['Warning!\nUpdating will permenently remove your ', ...
        'Chebfun root directory,\n%s/, and all contents.\n'], installDir);
    yesno = '';
    while ( ~any(strcmp(yesno, {'yes', 'no'}))  )
        fprintf(2, 'Proceed (yes / no)? ');
        yesno = input('', 's');
        if ( strcmpi(yesno, 'no') )
            cd(startDir)
            return
        end
    end
    
    % Check again!
    yesno = '';
    while ( ~any(strcmp(yesno, {'yes', 'no'}))  )
        fprintf(2, 'Are you sure (yes / no)? ');
        yesno = input('', 's');
        if ( strcmpi(yesno, 'no') )
            cd(startDir)
            return
        end
    end
    
    % Check for a Git repo:
    if ( exist('.git', 'dir') )
        fprintf(2, ['Warning: Your Chebfun root directory appears to be ', ...
            'part of a Git repo.\n']);
        yesno = '';
        while ( ~any(strcmp(yesno, {'yes', 'no'}))  )
            fprintf(2, 'Do you still want to continue (yes / no)? ');
            yesno = input('', 's');
            if ( strcmpi(yesno, 'no') )
                cd(startDir)
                return
            end
        end
    end
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try % Wrap everything in a try-catch statement to avoid a disaster.
        
    dstr = datestr(now, 'yyyymmddHHMMSS');
    if ( backup ) %#ok<*UNRCH>
        % Make a back up, so that we can recover things if needs be.
        filename = ['chebfun_backup_', dstr, '.zip']; 
        fprintf('Saving backup .zip to %s/ ...\n', fullfile(startDir, filename));
        warnState = warning('off', 'MATLAB:zip:archiveName');
        zip(filename, '*', pwd);
        warning(warnState);
        movefile(filename, startDir, 'f')
    end
    
    % Get the requested version of Chebfun:
    if ( strcmp(version, 'release') )
        zipDirName = 'chebfun-master';
        disp('Downloading .zip file of latest stable release ...')
        unzip([gitHubURL, 'archive/master.zip'], installDir);
    else
        zipDirName = 'chebfun-development';
        disp('Downloading .zip file of latest development version ...')
        unzip([gitHubURL, 'archive/development.zip'], installDir);
    end 
    
    if ( exist(fullfile(installDir, 'tmp'), 'dir') )
        rmdir(fullfile(installDir, 'tmp'), 's');
    end
        
    disp('Removing old files ...')   
    d = dir;
    mkdir(fullfile(installDir, 'tmp'));
    for k = 1:numel(d)
        % Remove subdirectories:
        if ( any(strcmp(d(k).name, {zipDirName, '.', '..'})) )
            % Don't remove downloaded zip, or current directory.
            continue
        end
        if ( d(k).isdir )
            rmdir(fullfile(installDir, d(k).name), 's');
        else
            % Move files to a tmp dir:
            movefile(fullfile(installDir, d(k).name), ...
                fullfile(installDir, 'tmp'), 'f');
        end
    end
    % remove the tmp dir:
    rmdir(fullfile(installDir, 'tmp'), 's');

    disp('Extracting .zip file ...')
    movefile(fullfile(installDir, zipDirName, '*'), installDir, 'f')
    rmdir(fullfile(installDir, zipDirName), 's')
    % Save the path:
    disp('Saving path ...')
    addpath(installDir)
    if ( savepath() )
        warning('CHEBFUN:TRUNKDIR:updateChebfun:savepath', ...
            'Unable to save path. Check with with your local sysadmin.');
    end
    
    % Return to the starting directory:
    cd(startDir)

    % Tidy up:
    close all
    clear all
    clear classes
    
    % Display outgoing info:
    disp('Installation complete!')
    
%     % Suggest they run CHEBTEST():
%     if ( usejava('jvm') && usejava('desktop') )
%         str = '<a href="matlab: chebtest">chebtest</a>';
%     else
%         str = 'chebtest';
%     end
%     disp(['We recommend you run ' str ' to ensure everything is working.'])
    
catch ME
       
    Return to the starting directory:
    cd(startDir)
    
    Something went wrong. Throw an error:
    disp('Update failed.');
    rethrow(ME)
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
