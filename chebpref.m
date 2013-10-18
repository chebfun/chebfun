classdef chebpref
% TODO:  Add documentation for this file.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    properties
        maxLength
        enableBreakpointDetection
        breakpointPrefs
        domain
        enableSingularityDetection
        singPrefs
        tech
        techPrefs
    end

    methods

        function p = chebpref(q)
            if ( (nargin == 1) && isa(q, 'chebpref') )
                p = q;
                return
            elseif ( nargin < 1 )
                q = struct();
            end

            % Initialize default preference values.
            p.maxLength = 65536;
            p.enableBreakpointDetection = false;
                p.breakpointPrefs.splitLength = 128;
                p.breakpointPrefs.splitMaxLength = 6000;
            p.domain = [-1 1];
            p.enableSingularityDetection = false;
                p.singPrefs.exponentTol = 1.1*1e-11;
                p.singPrefs.maxPoleOrder = 20;
            p.tech = 'chebtech';
            p.techPrefs = struct();

            % Copy fields from q, placing unknown ones in techPrefs and merging
            % incomplete substructures.
            for (field = fieldnames(q).')
                if ( isprop(p, field{1}) )
                    if ( isstruct(p.(field{1})) )
                        p.(field{1}) = chebpref.mergePrefs(p.(field{1}), ...
                            q.(field{1}));
                    else
                        p.(field{1}) = q.(field{1});
                    end
                else
                    p.techPrefs.(field{1}) = q.(field{1});
                end
            end
        end

        function out = subsref(p, ind)
        %SUBSREF   Subscripted referencing for CHEBPREF.
        %   P.PROP, where P is a CHEBPREF object, returns the value of the
        %   CHEBPREF property PROP stored in P.  If PROP is not a CHEBPREF
        %   property, P.TECHPREFS.PROP will be returned instead.  If PROP is
        %   neither a CHEBPREF property nor a field in P.TECHPREFS, an error
        %   will be thrown.
        %
        %   For access to fields PROP of TECHPREFS that have the same name as a
        %   CHEBPREF property, use the syntax P.TECHPREFS.PROP.
        %
        %   CHEBPREF does not support any other subscripted referencing types,
        %   including '()' and '{}'.
            switch ( ind(1).type )
                case '.'
                    if ( isprop(p, ind(1).subs) )
                        out = p.(ind(1).subs);
                    else
                        out = p.techPrefs.(ind(1).subs);
                    end

                    if ( numel(ind) > 1 )
                        out = subsref(out, ind(2:end));
                    end
                otherwise
                    error('CHEBTECH:subsref:badType', ...
                        'Invalid subscripted reference type.')
            end
        end

        function p = subsasgn(p, ind, val)
        %SUBSASGN   Subscripted assignment for CHEBPREF.
        %   P.PROP = VAL, where P is a CHEBPREF object, assigns the value VAL
        %   to the CHEBPREF property PROP stored in P.  If PROP is not a
        %   CHEBPREF property, the assignment will be made to P.TECHPREFS.PROP
        %   instead.
        %
        %   To assign to fields PROP of TECHPREFS that have the same name as a
        %   CHEBPREF property, use the syntax P.TECHPREFS.PROP = VAL.
        %
        %   CHEBPREF does not support any other subscripted assignment types,
        %   including '()' and '{}'.
            switch ( ind(1).type )
                case '.'
                    if ( isprop(p, ind(1).subs) )
                        p = builtin('subsasgn', p, ind, val);
                    else
                        p.techPrefs = builtin('subsasgn', p.techPrefs, ind, ...
                            val);
                    end
                otherwise
                    error('CHEBTECH:subsasgn:badType', ...
                        'Invalid subscripted assignment type.')
            end
        end

    end

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
