classdef chebpref
% TODO:  Add documentation for this file.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    methods ( Static = true )

        function p = mergePrefs(p, q, map)
        %MERGEPREFS   Merge preference structures.
        %   P = MERGEPREFS(P, Q) "merges" the preference structure Q into the
        %   preference structure P by replacing the contents of fields in P
        %   with those of identically-named fields in Q.  If Q has a field
        %   whose name does not match any of those in P, it is added to P.
        %
        %   P = MERGEPREFS(P, Q, MAP) does the same but uses the structure MAP
        %   to "translate" the names of fields of Q into names of fields of P.
        %   If Q has a field FIELD and MAP has a field of the same name, then
        %   the value of P.(MAP.FIELD) will be replaced by the contents of
        %   Q.FIELD.  If P does not have a field matching the string stored in
        %   MAP.FIELD, one will be added to P.
        %
        %   Note that this function operates on plain MATLAB structures, not
        %   CHEBPREF objects.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Developer notes:
        %  - This function is a helper function intended for use by "technology"
        %    objects (usually subclasses of SMOOTHFUN) for managing their
        %    preferences.  See CHEBTECH.PREF for an illustration.
        %  - The second syntax is useful, e.g., if Q contains abstractly-named
        %    preferences which may have a better name within the specific
        %    context of the tech object whose preferences are stored in P.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ( nargin < 3 )
                map = struct();
            end

            for (field = fieldnames(q).')
                if ( isfield(map, field{1}) )
                    p.(map.(field{1})) = q.(field{1});
                else
                    p.(field{1}) = q.(field{1});
                end
            end
        end

    end
end
