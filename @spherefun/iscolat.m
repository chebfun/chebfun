function b = iscolat( f )
% B = ISCOLAT( d ) Determines if the latitudinal variable of the SPHEREFUN
% is measured in co-latitude.
%
% Returns true of the domain of the latitudinal variable of the SPHEREFUN
% is measured in co-latitude, i.e. the North pole = 0 radians latitude,
% South Pole = pi radians latitude.

if isa( f, 'spherefun' )
    domain = f.domain;
    domain = domain(3:4);
elseif isnumeric(f) && numel(f) == 2
    domain = f;
else
    error('SPHEREFUN:iscolat:unknown',['Input argument of type %s unknown. ' ...
    'Input should be a spherefun or a vector with 2 elements, of which the ' ...
    'entries are the domain for latitude'],class(f));
end

% Check latitudinal coordinate system to set the elevation angle
% appropriately.
colat = [0 pi]; % Colatitude (doubled up)
b = all( (domain-colat) == 0 );

end

