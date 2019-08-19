function data = dispData(f)
%DISPDATA   Useful information for DISPLAY at higher levels.
%   DATA = DISPDATA(F) extracts useful information from the given CLASSICFUN 
%   F and the information will be used by DISPLAY at higher levels. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

data = dispData(f.onefun);

if ( isa(f.mapping, 'singularMapping') )
    data{1}(end+1).name = '  singmap params ';
    params = f.mapping.params;
    data{1}(end).data = ['  [' num2str(params(1), '%2.2g,') '     ' ...
            num2str(params(2), '%2.2g') ']  '];
end

% More information for F can be appended to DATA:

end
