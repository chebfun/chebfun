function val = get( f, propName )
%GET       GET method for DISKFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the DISKFUN object F. Valid entries for the string PROP are:
%    'domain'
%    'cols'
%    'rows' 
%    'pivotValues'
%    'pivotLocations'
%    'pivotIndices'
%    'idxPlus'
%    'idxMinus'
%    'nonZeroPoles'

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Loop through an array of DISKFUN objects.
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
    case 'idxPlus'
        val = f.idxPlus;
    case 'idxMinus'
        val = f.idxMinus;
    case 'nonZeroPoles'
        val = f.nonZeroPoles;
    otherwise
        error('CHEBFUN:DISKFUN:get:propName', ...
            [propName,' is not a valid DISKFUN property.'])
end
