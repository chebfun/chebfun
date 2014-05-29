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
        
        function printDescription(fid, expInfo)
            
            % Extract info from the expInfo struct:
            deInput = expInfo.deInput;
            indVarNameSpace = expInfo.indVarNameSpace;
            dom = expInfo.dom;
            bcInput = expInfo.bcInput;
            periodic = expInfo.periodic;
            
            % Print a description of the BVP:
            fprintf(fid, '%% Solving\n');
            for k = 1:numel(deInput)
                fprintf(fid, '%%   %s,\n', deInput{k});
            end
            fprintf(fid, '%% for %s in %s', indVarNameSpace, dom);
            if ( ~isempty(bcInput{1}) )
                fprintf(fid, ', subject to\n%%');
                for k = 1:numel(bcInput)
                    fprintf(fid, '   %s', bcInput{k});
                    if ( (k ~= numel(bcInput)) && (numel(bcInput) > 1) )
                        fprintf(fid, ',\n%%');
                    end
                end
                fprintf(fid, '.\n');
            elseif periodic
                fprintf(fid, ', subject to periodic boundary conditions.\n\n');
            else
                fprintf(fid, '.\n');
            end
            
        end
        
        printSetup(fid, expInfo, guifile)
    
    end
end