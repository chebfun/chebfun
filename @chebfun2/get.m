function val = get( f, propName )
%GET   GET method for CHEBFUN2 class.
%
% P = GET(F, PROP) returns the property P specified in the string PROP from
% the chebfun F. Valid entries for the string PROP are:
%   'DOMAIN'
%   'COLS'
%   'ROWS' 
%   'PIVOTVALUES'
%   'PIVOTLOCATIONS'

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Loop through an array of chebfun2 objects.
if ( numel(f) > 1 )
    val = cell(numel(f));
    for k = 1:numel(f)
        val{k} = get(f(k), propName);
    end
    return
end

%%
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
    otherwise
        error('CHEBFUN2:get:propnam',[propName,' is not a valid chebfun2 property.'])
end