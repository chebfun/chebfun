classdef chebguiExporter
%CHEBGUIEXPORTER   Export a problem from CHEBGUI.
%   This is a an abstract class that defines interfaces for export problems from
%   CHEBGUI to .m-files, to the workspace, or to a .chebgui file. It is not
%   intended to be called directly by the end user.
%
%   See also CHEBGUI, CHEBGUIBVPEXPORTER, CHEBGUIEIGEXPORTER,
%   CHEBGUIPDEEXPOERTER.

% Developers note:
%   The CHEBGUICONTROLLER class defines a number of abstract methods, used to
%   export problems from CHEBGUI. In v4, this functionality used to live in the
%   @chebgui folder, but to increase modularity, it has been spun off to its own
%   class.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    properties ( Access = public )
       
        % Default file name to be saved to
        defaultFileName
        
    end

    methods ( Abstract = true )
        
        % Export to .m file
        chebgui2mfile(exporter, guifile, pathname, filename)
        
    end
    
    
    methods( Abstract = true, Static = true )
        % Make a CHEBGUIEXPORTER. (Constructor shortcut)
        e = make(varargin);
        
    end
    
    methods ( Static = true )
        function obj = constructor(type)
            % Constructor for the CHEBGUIEXPORTER class.
            
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
    end
    
    methods ( Access = public )
        
        % Export to an .m-file.
        toFile(e, guifile)
        
    end
    
end