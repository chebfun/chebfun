function prefs = pref(varargin)
%PREF   Preference settings for SMOOTHFUN.
%   SMOOTHFUN.PREF(PREFNAME) returns the value corresponding to the preference
%   named in the string PREFNAME.
%
%   P = SMOOTHFUN.PREF returns a structure P with a field P.SMOOTHFUN which contains
%   the default SMOOTHFUN preferences as fields/values. This structure may be
%   passed to the SMOOTHFUN constructor.
%
%   P = SMOOTHFUN.PREF(P) will check to see whether the input preference structure
%   P already has a SMOOTHFUN field. If it does not, one is appended.
%
%   P = SMOOTHFUN.PREF('PREFNAME1', VAL1, 'PREFNAME2', VAL2, ...) returns the same
%   structure as above, but with the default SMOOTHFUN preferences 'PREFNAME1',
%   'PREFNAME2', replaced by the values in VAL1, VAL2, etc.
%
%   P = SMOOTHFUN.PREF(P, 'PREFNAME1', VAL1, 'PREFNAME2', VAL2, ...) appends a
%   SMOOTHFUN preference field to P if required, and modifies the SMOOTHFUN
%   properties 'PREFNAME1' and 'PREFNAME2'.
%
%   Note that no checking of either the input PREFNAMEs or VALs takes place.
%
%   SMOOTHFUN PREFERENCES (case sensitive)
%
%     eps          -  Relative tolerance used in construction and subsequent
%      [2^-52]        operations.

% See also CHEBTECH, CHEBTECH1, CHEBTECH2

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

classname = 'smoothfun';

% Has a preference structure been passed?
if ( numel(varargin) > 0 && isstruct(varargin{1}) ) % Yes
    prefs = varargin{1};
    varargin(1) = [];
else                                                % No
    prefs = struct();
end

% If it was, did it have a field for this class?
if ( isfield(prefs, classname) )  % It does, so either:
    if ( numel(varargin) == 0 )
        return                    % a) No props to change, return
    else
        p = prefs.(classname);    % b) Grab prefs for this class.
    end
else                              % No prefs found for this class, so make some:
    p.eps    = 2^-52;
end
% p is now the preference substructure relating to the current class.

if ( isfield(prefs,'misc') ) 
    q = prefs.misc;
else
    q = struct();
end
% q is now the preference substructure for MISC preferences.

for names = fieldnames(q)'
    if isfield(p,names)
        field = names{:};
        p.(field) = q.(field);
        q = rmfield(q,field);
    end
end

% Two preference structures were passed. Copy matching fieldnames:
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

% A single property has been queried, so return this property:
if ( numel(varargin) == 1 )
    prefs = p.(varargin{1});
    return
end

% Property names have been passed, so alter/add <current class>/MISC properties.
for k = 1:2:numel(varargin)
    if ( isfield(p, varargin{k}) )
        p.(varargin{k}) = varargin{k+1};
    else
        q.(varargin{k}) = varargin{k+1};
    end
end

% Append current class preferences to the preference structure prefs for output:
prefs.(classname) = p;

% Append MISC preferences to the preference structure prefs for output:
prefs.misc = q;

end