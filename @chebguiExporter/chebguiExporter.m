classdef chebguiExporter
%CHEBGUIEXPORTER   Class to export a problem from CHEBGUI.
%   This is a an abstract class that defines interfaces for export problems from
%   CHEBGUI to .m-files, to the workspace, or to a .chebgui file. It is not
%   intended to be called directly by the end user.
%
% See also CHEBGUI, CHEBGUIBVPEXPORTER, CHEBGUIEIGEXPORTER, CHEBGUIPDEEXPOERTER.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%   The CHEBGUIEXPORTER class defines a number of abstract methods, used to
%   export problems from CHEBGUI. In v4, this functionality used to live in the
%   @chebgui folder, but to increase modularity, it has been spun off to its own
%   class and subclasses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Abstract = true )
        
        % Default file name to be saved to
        defaultFileName
        
        % Description of exported .m files
        description
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Abstract = true, Static = true )
        
        % Extract information from the CHEBGUI object to a struct:
        expInfo = exportInfo(guifile)
        
        % Print problem description:
        printDescription(fid, expInfo)
        
        % Print options for solving the problem:
        printOptions(fid, expInfo)
        
        % Print lines for setting up the problem:
        printSetup(fid, expInfo, guifile)
        
        % Print the actual lines for calling the solver method:
        printSolver(fid, expInfo)
        
        % Print steps taken after the solver finishes:
        printPostSolver(fid, expInfo)
        
        % Export to a .mat-file
        toMat(handles)
        
        % Export to workspace
        toWorkspace(handles)
        
        % Export to workspace, solution only
        toWorkspaceSolutionOnly(handles)

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods ( Access = public, Static = false )
        
        function toFile(exporter, guifile, fileName, pathName)
        %TOFILE   Export a CHEBGUI to an m-file
        %   TOFILE(EXPORTER, GUIFILE, FILENAME, PATHNAME) exports the CHEBGUI
        %   object GUIFILE to the file 'PATHNAME/FILENAME' using the
        %   CHEBGUIEXPORTER object EXPORTER (of BVP, EIG, or PDE kind).

        % Concatenate the pathName and the fileName to get the full path:
        fullFileName = fullfile(pathName, fileName);

        try
            % Open a stream to write to a file:
            fid = fopen(fullFileName, 'wt');

            % Extract the necessary info for export to an .m file from the
            % GUIFILE object:
            expInfo = exporter.exportInfo(guifile);

            % Write the header information:
            writeHeader(exporter, fid, fileName)

            % Print description of the problem:
            exporter.printDescription(fid, expInfo)

            % Print the problem set-up:
            exporter.printSetup(fid, expInfo, guifile)

            % Print the options set-up:
            exporter.printOptions(fid, expInfo)

            % Print the solution step:
            exporter.printSolver(fid, expInfo)

            % Print the post-solution process:
            exporter.printPostSolver(fid, expInfo)

            % Close the file writing stream:
            fclose(fid);
        catch ME
            % Make sure to tidy up first
            fclose(fid);

            % Rethrow the error:
            rethrow(ME)
        end
        
        end
        
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        function obj = constructor(type)
        %CONSTRUCTOR   Constructor for the CHEBGUIEXPORTER class.
        %   OBJ = CONSTRUCTOR(TYPE), where TYPE = 'bvp', 'eig' or 'pde' returns
        %   an object that is a concrete implementation of CHEBGUIEXPORTER
        %   object OBJ. The type of OBJ is determined by the TYPE input
        %   argument.
        %
        %   This method allows one to construct concrete subclasses of the
        %   CHEBGUIEXPORTER class.
            
            % Return the desired type of CHEBGUIEXPORTER
            switch type
                case 'bvp'
                    obj = chebguiExporterBVP;
                    
                case 'eig'
                    obj = chebguiExporterEIG;

                case 'ivp'
                    obj = chebguiExporterIVP;
                    
                case 'pde'
                    obj = chebguiExporterPDE;
                    
                otherwise
                    error('CHEBFUN:CHEBGUIEXPORTER:chebguiExporter:constructor', ...
                        'Unknown type for CHEBGUIEXPORTER constructor.')
            end
            
        end
        
        function toChebgui(guifile)
        %TOCHEBGUI   Export CHEBGUI to a CHEBGUI object in the workspace.
            
            % Setup for input dialog
            prompt = 'Enter the name of the chebgui variable:';
            name = 'Export GUI';
            numlines = 1;
            defaultAnswer ='chebg';
            options.Resize ='on';
            options.WindowStyle ='modal';
            
            % Get user input
            answer = inputdlg(prompt, name, numlines, {defaultAnswer}, options);
            
            % If user did not press cancel, assign to a variable in the base
            % workspace.
            if ( ~isempty(answer) )
                assignin('base', answer{1}, guifile);
            end
        end

        
        function str = prettyPrintFevalString(str, varnames)
        %PRETTYPRINTFEVALSTRING   Tidy up strings printed to .m files.
            
            % Do various string manipulations to make the output look nicer.
            for k = 1:numel(varnames)
                oldstr = ['feval(' varnames{k} ','];
                newstr = [varnames{k} '('];
                str = strrep(str, oldstr, newstr);
                oldstr = [varnames{k} '(''end'''];
                newstr = [varnames{k} '(end'];
                str = strrep(str, oldstr, newstr);
                oldstr = [varnames{k} '(''right'''];
                newstr = [varnames{k} '(end'];
                str = strrep(str, oldstr, newstr);
                oldstr = [varnames{k} '(''start'''];
                newstr = [varnames{k} '(' varnames{k} '.ends(1)'];
                str = strrep(str, oldstr, newstr);
                oldstr = [varnames{k} '(''left'''];
                newstr = [varnames{k} '(' varnames{k} '.ends(1)'];
                str = strrep(str, oldstr, newstr);
            end
            
        end
        
        function disc = discOption(isPeriodic, dom, opt)
            % Do we want to use TRIGCOLLOC for the discretization? This will be
            % the case if we're solving a periodic problem, and have
            % 'collocation' specified as the discretization option. However, we
            % don't want to use TRIGCOLLOC if there are breakpoints in the
            % domain.
            if ( isPeriodic && strcmp(opt, 'collocation') &&  ...
                    ( length(str2num(dom)) == 2 ) )
                disc = 'periodic';
            else
                disc = opt;
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE PROTECTED METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = protected )
        
        function writeHeader(e, fid, filename)
            %WRITEHEADER   The first few lines of exported .m-files.
            %
            % Inputs:
            %   E:          a CHEBGUIEXPORT object
            %   FID:        a valid file ID
            %   FILENAME:   name of the file being written to
            
            % Access usernames in different ways, depending on whether are
            % running on PCs or not:
            if ( ispc )
                userName = getenv('UserName');
            else
                userName = getenv('USER');
            end
            
            % Print first few lines of the .m-file:
            fprintf(fid, ['%%%% %s -- an executable m-file for solving ', ...
                '%s\n'], filename, e.description);
            fprintf(fid, ['%% Automatically created in CHEBGUI ', ...
                'by user %s.\n'], userName);
            fprintf(fid, '%% Created on %s', ...
                datestr(now, 'mmmm dd, yyyy at HH:MM.\n\n'));
            
            fprintf(fid, '%%%% Problem description.\n');
            
        end
        
    end
    
end
