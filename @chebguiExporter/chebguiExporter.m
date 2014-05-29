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
        
    methods ( Abstract = true)
        
        % Export to .m file
        chebgui2mfile(exporter, guifile, pathname, filename)
        
    end
    
    
    methods( Abstract = true, Static = true )
        % Make a CHEBGUIEXPORTER. (Constructor shortcut)
        e = make(varargin);
        
    end
end