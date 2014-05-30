classdef chebguiExporterBVP < chebguiExporter
    %CHEBGUIEXPORTERBVP   Export a BVP from CHEBGUI.
    %   This is a an concrete implementation of the class CHEBGUIEXPORTER, which
    %   exports BVPs from CHEBGUI to .m-files, to the workspace, or to a .chebgui
    %   file. It is not intended to be called directly by the end user.
    %
    %   See also CHEBGUI, CHEBGUIEXPORTER.
    
    % Developers note:
    %   The CHEBGUICONTROLLER class defines a number of abstract methods, used to
    %   export problems from CHEBGUI. In v4, this functionality used to live in the
    %   @chebgui folder, but to increase modularity, it has been spun off to its own
    %   class.
    
    % Copyright 2014 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    methods (Access = public)
        
        function A = chebguiExporterBVP(varargin)
            % Set default file name:
            A.defaultFileName = 'bvp.m';
            
            % Description for printing to .m files:
            A.description = 'a boundary-value problem';
        end
        
    end
    
    
    methods ( Static = true )
        
        function e = make(varargin)
            % Factory method.
            e = chebguiExporterBVP(varargin{:});
        end
        
        function printSolver(fid, expInfo)
            %PRINTSOLVER    Print the solution step when exporting
            fprintf(fid,'\n%%%% Solve the problem!');
            fprintf(fid, ['\n%% Here, we call the solvebvp() method ' ...
                '(which offers the same functionality \n%% as nonlinear '...
                'backslash, but with more customizability).\n']);
            fprintf(fid, 'u = solvebvp(N, rhs, options);\n');
        end
        
        function printPostSolver(fid, expInfo)
            %PRINTPOSTSOLVER    Print commands after solution has been found
            
            % Extract information from the EXPINFO struct
            allVarNames = expInfo.allVarNames;
            indVarNameSpace = expInfo.indVarNameSpace;
            
            % Print commands that will create a plot of the solution obtained:
            fprintf(fid, '\n%%%% Create a plot of the solution.\n');
            
            fprintf(fid, ['figure\nplot(u,''LineWidth'',2)\n', ...
                'title(''Final solution''), xlabel(''%s'')'], indVarNameSpace);
            if ( numel(allVarNames) == 1 )
                % Scalar problem:
                fprintf(fid, ', ylabel(''%s'')', allVarNames{:});
            else
                % Coupled system. Create a legend.
                leg = '';
                for k = 1:numel(allVarNames)-1
                    leg = [leg '''' allVarNames{k} '''' ','];
                end
                leg = [leg '''' allVarNames{k+1} ''''];
                fprintf(fid, ', legend(%s)\n', leg);
            end
            
            
        end
        
        printDescription(fid, expInfo)
        
        printSetup(fid, expInfo, guifile)
        
        printOptions(fid, expInfo, guifile)
    end
end