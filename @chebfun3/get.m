function val = get(f, propName)
%GET   GET method for CHEBFUN3 class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the CHEBFUN F. Valid entries for the string PROP are:
%    'COLS'
%    'ROWS' 
%    'TUBES'
%    'CORE'
%    'DOMAIN'

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