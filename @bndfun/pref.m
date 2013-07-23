function prefs = pref(varargin)
% [TODO]: Do we actually need this preference method? If so, it needs
% documentation.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

classname = 'bndfun';

% Has a preference structure been passed?
if ( numel(varargin) > 0 && isstruct(varargin{1}) ) % Yes
    prefs = varargin{1};
    varargin(1) = [];
else                                                % No
    prefs = struct();
end

% If it was, did it have a BNDFUN field?
if ( isfield(prefs, classname) )  % It does, so either:
    if ( numel(varargin) == 0 )
        return                    % a) No props to change, return
    else
        p = prefs.(classname);    % b) Grab bndfun prefs
    end
else                              % No bndfun prefs found, so make some:
    p.eps         = 2^-52;
    p.domain      = [-1, 1];
end
% p is now the preference substructure relating to the FUN class.

if ( isfield(prefs,'misc') ) 
    q = prefs.misc;
else
    q = struct();
end
% q is now the preference substructure for MISC preferences

for names = fieldnames(q)'
    if isfield(p,names)
        field = names{:};
        p.(field) = q.(field);
        q = rmfield(q,field);
    end
end

% Two preference structures were passed. Copy matching fieldnames.
if ( numel(varargin) > 0 && isstruct(varargin{1}) ) % Yes
    r = varargin{1};
    varargin(1) = [];
    for names = fieldnames(r).'
        if isfield(p,names)
            field = names{:};
            p.(field) = r.(field);
        end
    end
end

% A single property has been queried, so return this property
if ( numel(varargin) == 1 )
    prefs = p.(varargin{1});
    return
end

% Property names have been passed, so alter/add BNDFUN/MISC properties.
for k = 1:2:numel(varargin)
    if ( isfield(p, varargin{k}) )
        p.(varargin{k}) = varargin{k+1};
    else
        q.(varargin{k}) = varargin{k+1};
    end
end

% Append BNDFUN preferences to the preference structure prefs for output.
prefs.(classname) = p;

% Append MISC preferences to the preference structure prefs for output.
prefs.misc = q;

