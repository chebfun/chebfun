function val = get(f, propName)
%GET    GET method for SPHEREFUN class.
%   P = GET(F, PROPNAME) returns the property P specified in the string 
%   PROPNAME from the SPHEREFUN object F. Valid entries for the string 
%   PROPNAME are:
%    'DOMAIN'
%    'COLS'
%    'ROWS' 
%    'PIVOTVALUES'
%    'PIVOTLOCATIONS'
%    'BLOCKDIAG'
%    'IDXPLUS'
%    'IDXMINUS'
%    'NONZEROPOLES'

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Loop through an array of SPHEREFUN objects.
if ( numel(f) > 1 )
    val = cell(numel(f));
    for k = 1:numel(f)
        val{k} = get(f(k), propName);
    end
    return
end

% Get the properties.
switch ( propName )
    case 'domain'
        val = f.domain;
    case 'cols'
        val = f.cols;
    case 'rows'
        val = f.rows;
    case 'pivotValues'
        val = f.pivotValues;
    case 'pivotLocations'
        val = f.pivotLocations;
    case 'blockDiag'
        val = f.blockDiag;
    case 'idxPlus'
        val = f.idxPlus;
    case 'idxMinus'
        val = f.idxMinus;
    case 'nonZeroPoles'
        val = f.nonZeroPoles;
    otherwise
        error('CHEBFUN:SPHEREFUN:get:propName', ...
            [propName,' is not a valid SPHEREFUN property.'])
end