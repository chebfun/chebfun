function export(guifile, handles, exportType)

% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

problemType = guifile.type;

% handles and exportType are empty if user calls export from command window.
% Always export to .m file in that case.
if ( nargin == 1 )
    exportType = '.m';
end

% Exporting a BVP
if ( strcmp(problemType, 'bvp') )
    switch exportType
        case 'Workspace'
            
            exporter = chebguiExporterBVP;
            exporter.toWorkspace(handles)
            
        case 'Cancel'
            return;
    end
    
% Exporting a PDE
elseif ( strcmp(problemType,'pde') )
    switch exportType
        case 'Workspace'
            
            exporter = chebguiExporterPDE;
            exporter.toWorkspace(handles)
            
            
        case 'Cancel'
            return;
    end

% Exporting an EIG problem.
else
    switch exportType
        case 'Workspace'
            
            exporter = chebguiExporterEIG;
            exporter.toWorkspace(handles)
            
           
        case 'Cancel'
            return;
    end
end
