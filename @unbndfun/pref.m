function prefs = pref(varargin)

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

classname = 'unbndfun';

% Has a preference structure been passed?
if ( numel(varargin) > 0 && isstruct(varargin{1}) ) % Yes
    prefs = varargin{1};
    varargin(1) = [];
else                                                % No
    prefs = struct();
end

% If it was, has it have a unbndfun field already?
if ( isfield(prefs, classname) )  % It does, so either:
    if ( numel(varargin) == 0 )
        return                    % a) No props to change, return
    else
        p = prefs.(classname);    % b) Grab unbndfun prefs
    end
else                              % No unbndfun prefs found, so make some:
    p.eps         = 2^-52;
    p.domain      = [-1, 1];     %[TODO]: default domain
end
% p is now the preference substructure relating to the UNBNDFUN class.

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

% Property names have been passed, so alter/add UNBNDFUN/MISC properties.
for k = 1:2:numel(varargin)
    if ( isfield(p, varargin{k}) )
        p.(varargin{k}) = varargin{k+1};
    else
        q.(varargin{k}) = varargin{k+1};
    end
end

% Append UNBNDFUN preferences to the prefence structure prefs for output.
prefs.(classname) = p;

% Append MISC preferences to the prefence structure prefs for output.
prefs.misc = q;