classdef chebguiExporterEIG < chebguiExporter
    %CHEBGUIEXPORTEREIG   Export an EIG problem from CHEBGUI.
    %   This is a an concrete implementation of the class CHEBGUIEXPORTER, which
    %   exports EIG problems from CHEBGUI to .m-files, to the workspace, or to a
    %   .chebgui file. It is not intended to be called directly by the end user.
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
        
        function A = chebguiExporterEIG(varargin)
            % Set default file name:
            A.defaultFileName = 'bvpeig.m';
            
            % Description for printing to .m files:
            A.description = 'an eigenvalue problem';
        end
        
    end
    
    
    methods ( Static = true )
        
        function e = make(varargin)
            % Factory method.
            e = chebguiExporterEIG(varargin{:});
        end
        
        function printOptions(fid, expInfo)
            K = expInfo.K;
            
            fprintf(fid, '\n%% Number of eigenvalue and eigenmodes to compute.\n');
            fprintf(fid, 'k = %s;\n', K);
            
            fprintf(fid, '\n%% Number of eigenvalue and eigenmodes to compute.\n');
            fprintf(fid, 'k = %s;\n', K);
            
        end
        
        function printSolver(fid, expInfo)
            sigma = expInfo.sigma;
            generalized = expInfo.generalized;
            
            fprintf(fid, '\n%%%% Solve the eigenvalue problem.\n');
            if ( ~generalized )
                if ( ~isempty(sigma) )
                    fprintf(fid, '[V, D] = eigs(N, k, %s);\n', sigma);
                else
                    fprintf(fid, '[V, D] = eigs(N, k);\n');
                end
            else
                if ( ~isempty(sigma) )
                    fprintf(fid, '[V, D] = eigs(N, B, k, %s);\n', sigma);
                else
                    fprintf(fid, '[V, D] = eigs(N, B, k);\n');
                end
            end
        end
        
        function printPostSolver(fid, expInfo)
            
            allVarNames = expInfo.allVarNames;
            indVarName = expInfo.indVarName;
            allVarString = expInfo.allVarString;
            
            fprintf(fid, '\n%%%% Plot the eigenvalues.\n');
            fprintf(fid, 'D = diag(D);\n');
            fprintf(fid, 'figure\n');
            fprintf(fid, 'plot(real(D), imag(D), ''.'', ''markersize'', 25)\n');
            fprintf(fid, 'title(''Eigenvalues''); xlabel(''real''); ylabel(''imag'');\n');
            
            if ( ischar(allVarNames) || (numel(allVarNames) == 1) )
                fprintf(fid, '\n%% Plot the eigenmodes.\n');
                fprintf(fid, 'figure\n');
                fprintf(fid, 'plot(real(V), ''linewidth'', 2);\n');
                fprintf(fid, 'title(''Eigenmodes''); xlabel(''%s''); ylabel(''%s'');', ...
                    indVarName{1}, allVarString);
            end
        end
        
        printDescription(fid, expInfo)
        
        printSetup(fid, expInfo, guifile)
        

    end
end