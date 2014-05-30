classdef chebguiExporterPDE < chebguiExporter
    %CHEBGUIEXPORTERPDE   Export a PDE from CHEBGUI.
    %   This is a an concrete implementation of the class CHEBGUIEXPORTER, which
    %   exports PDEs from CHEBGUI to .m-files, to the workspace, or to a .chebgui
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
        
        function A = chebguiExporterPDE(varargin)
            % Set default file name:
            A.defaultFileName = 'pde.m';
            
            % Description for printing to .m files:
            A.description = 'a partial differential equation.';
        end
        
    end
    
    
    methods ( Static = true )
        
        function e = make(varargin)
            % Factory method.
            e = chebguiExporterPDE(varargin{:});
        end
        
        printDescription(fid, expInfo)
        
        printOptions(fid, expInfo)
        
        printSetup(fid, expInfo, guifile)
        
        function printPostSolver(fid, expInfo)
            sol = expInfo.sol;
            indVarName = expInfo.indVarName;
            deInput = expInfo.deInput;
            s = expInfo.s;
            
            % plotting
            if ( numel(deInput) == 1 )
                fprintf(fid, '\n%% Create plot of the solution.\n');
                %     fprintf(fid,'surf(%s,t,''facecolor'',''interp'')\n',sol);
                fprintf(fid, 'waterfall(%s,%s,''simple'',''linewidth'',2)\n', sol, ...
                    indVarName{2});
            else
                fprintf(fid, '\n%% Create plots of the solutions.\n');
                M = numel(deInput);
                for k = 1:numel(deInput)
                    fprintf(fid, 'subplot(1,%d,%d)\n', M, k);
                    fprintf(fid, 'waterfall(%s,%s,''linewidth'',2)\n', s{k}, indVarName{2});
                    fprintf(fid, 'xlabel(''%s''), ylabel(''%s''), title(''%s'')\n', ...
                        indVarName{1},indVarName{2},s{k});
                end
            end
        end
        
        
        function printSolver(fid, expInfo)
            indVarName = expInfo.indVarName;
            sol = expInfo.sol;
            sol0 = expInfo.sol0;
            deInput = expInfo.deInput;
            s = expInfo.s;
            
            fprintf(fid, ['\n%% Solve the problem using pde15s.\n']);
            fprintf(fid, '[%s, %s] = pde15s(pdefun, %s, %s, bc, opts);\n', indVarName{2}, ...
                sol, indVarName{2}, sol0);
            
            % Conver sol to variable names
            if ( numel(deInput) > 1 )
                fprintf(fid, '\n%% Recover variable names.\n');
                for k = 1:numel(s)
                    fprintf(fid, '%s = %s{%d};\n', s{k}, sol, k);
                end
            end
        end
    end
end