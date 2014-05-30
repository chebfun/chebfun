classdef chebguiExporter
    %CHEBGUIEXPORTER   Export a problem from CHEBGUI.
    %   This is a an abstract class that defines interfaces for export problems from
    %   CHEBGUI to .m-files, to the workspace, or to a .chebgui file. It is not
    %   intended to be called directly by the end user.
    %
    %   See also CHEBGUI, CHEBGUIBVPEXPORTER, CHEBGUIEIGEXPORTER,
    %   CHEBGUIPDEEXPOERTER.
    
    % Developers note:
    %   The CHEBGUIEXPORTER class defines a number of abstract methods, used to
    %   export problems from CHEBGUI. In v4, this functionality used to live in
    %   the @chebgui folder, but to increase modularity, it has been spun off to
    %   its own class and subclasses.
    
    % Copyright 2014 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    properties ( Abstract )
        
        % Default file name to be saved to
        defaultFileName
        
        % Description of exported .m files
        description
        
    end
    
    methods ( Abstract = true )
        
        % Extract information from the CHEBGUI object to a struct
        expInfo = exportInfo(e, guifile)
        
    end
    
    methods ( Abstract = true, Static = true )
        
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
        
    end
    
    methods ( Static = true )
        function obj = constructor(type)
            %CONSTRUCTOR    Constructor for the CHEBGUIEXPORTER class.
            %
            % This method allows us to construct concrete subclasses of the
            % CHEBGUIEXPORTER class.
            
            % Return the desired type of CHEBGUIEXPORTER
            switch type
                case 'bvp'
                    obj = chebguiExporterBVP;
                    
                case 'eig'
                    obj = chebguiExporterEIG;
                    
                case 'pde'
                    obj = chebguiExporterPDE;
                    
                otherwise
                    error('CHEBFUN:CHEBGUIEXPORTER:constructor', ...
                        'Unknown type for CHEBGUIEXPORTER constructor.')
            end
            
        end
        
        function toChebgui(guifile)
            %TOCHEBGUI
            %
            % Export the CHEBGUI shown in the figure to a CHEBGUI variable in
            % the workspace
            
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
        
        function str = prettyPrintFevalString(str,varnames)
            %PRETTYPRINTFEVALSTRING     Tidy up strings printed to .m files
            
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
        
    end
    
    methods ( Access = public )
        
        % Export to an .m-file.
        toFile(e, guifile, fileName, pathName)
        
    end
    
    methods ( Access = protected )
        function writeHeader(e, fid, filename)
            %WRITEHEADER        The first few lines of exported .m-files.
            
            % Access usernames in different ways, depending on whether are
            % running on PCs or not:
            if ( ispc )
                userName = getenv('UserName');
            else
                userName = getenv('USER');
            end
            
            % Print first few lines of the .m-file:
            fprintf(fid, ['%% %s -- an executable .m-file for solving ', ...
                '%s.\n'], filename, e.description);
            fprintf(fid, ['%% Automatically created from chebfun/chebgui ', ...
                'by user %s\n'], userName);
            fprintf(fid, '%% at %s on %s.\n\n', datestr(rem(now, 1), 13), ...
                datestr(floor(now)));
            
        end
    end
end