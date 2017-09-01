function val = get(f, propName)
%GET   GET method for the CHEBFUN3 class.
%   P = GET(F, PROP) returns the property P of the CHEBFUN3 object F 
%   specified in the string PROP. Valid entries for the string PROP are:
%   'COLS'
%   'ROWS' 
%   'TUBES'
%   'CORE'
%   'DOMAIN'

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Loop through an array of CHEBFUN3S objects.
if ( numel(f) > 1 )
    val = cell(numel(f));
    for k = 1:numel(f)
        val{k} = get(f(k), propName);
    end
    return
end

% Get the properties.
switch ( propName )
    case 'cols'
        val = f.cols;
    case 'rows'
        val = f.rows;
    case 'tubes'
        val = f.tubes;        
    case 'core'
        val = f.core;
    case 'domain'
        val = f.domain;
    otherwise
        error('CHEBFUN:CHEBFUN3:get:propName', ...
            [propName,' is not a valid CHEBFUN3 property.'])
end