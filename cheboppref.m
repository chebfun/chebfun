classdef cheboppref
    
    % See above for documentation.
    properties
        domain
        discretisation
        discPrefs
        enableBreakpointDetection
        breakpointPrefs
        maxTotalLength
        scale
    end

    methods

        function outPref = cheboppref(inPref)
            if ( (nargin == 1) && isa(inPref, 'cheboppref') )
                outPref = inPref;
                return
            elseif ( nargin < 1 )
                inPref = struct();
            end

            % Initialize default preference values.
            outPref.maxTotalLength = 2048;
            outPref.enableBreakpointDetection = false;
                outPref.breakpointPrefs.splitMaxLength = 129;
                outPref.breakpointPrefs.splitMaxTotalLength = 2048;
            outPref.domain = [-1 1];
            outPref.discretisation = @colloc2;
            outPref.discPrefs = struct();
                outPref.discPrefs.eps = 2^(-52);
                outPref.discPrefs.exactLength = NaN;
            outPrefs.scale = NaN;    

            % Copy fields from q, placing unknown ones in discPrefs and merging
            % incomplete substructures.
            for field = fieldnames(inPref).'
                if ( isprop(outPref, field{1}) )
                    if ( isstruct(outPref.(field{1})) )
                        outPref.(field{1}) = ...
                            chebpref.mergePrefs(outPref.(field{1}), ...
                            inPref.(field{1}));
                    else
                        outPref.(field{1}) = inPref.(field{1});
                    end
                else
                    outPref.discPrefs.(field{1}) = inPref.(field{1});
                end
            end
        end

        function out = subsref(pref, ind)
        %SUBSREF   Subscripted referencing for CHEBPREF.
        %   P.PROP, where P is a CHEBPREF object, returns the value of the
        %   CHEBPREF property PROP stored in P.  If PROP is not a CHEBPREF
        %   property, P.DISCPREFS.PROP will be returned instead.  If PROP is
        %   neither a CHEBPREF property nor a field in P.DISCPREFS, an error
        %   will be thrown.
        %
        %   For access to fields PROP of DISCPREFS that have the same name as a
        %   CHEBPREF property, use the syntax P.DISCPREFS.PROP.
        %
        %   CHEBPREF does not support any other subscripted referencing types,
        %   including '()' and '{}'.
            switch ( ind(1).type )
                case '.'
                    if ( isprop(pref, ind(1).subs) )
                        out = pref.(ind(1).subs);
                    else
                        out = pref.discPrefs.(ind(1).subs);
                    end

                    if ( numel(ind) > 1 )
                        out = subsref(out, ind(2:end));
                    end
                otherwise
                    error('CHEBTECH:subsref:badType', ...
                        'Invalid subscripted reference type.')
            end
        end

        function pref = subsasgn(pref, ind, val)
        %SUBSASGN   Subscripted assignment for CHEBPREF.
        %   P.PROP = VAL, where P is a CHEBPREF object, assigns the value VAL
        %   to the CHEBPREF property PROP stored in P.  If PROP is not a
        %   CHEBPREF property, the assignment will be made to P.DISCPREFS.PROP
        %   instead.
        %
        %   To assign to fields PROP of DISCPREFS that have the same name as a
        %   CHEBPREF property, use the syntax P.DISCPREFS.PROP = VAL.
        %
        %   CHEBPREF does not support any other subscripted assignment types,
        %   including '()' and '{}'.
            switch ( ind(1).type )
                case '.'
                    if ( isprop(pref, ind(1).subs) )
                        pref = builtin('subsasgn', pref, ind, val);
                    else
                        pref.discPrefs = builtin('subsasgn', pref.discPrefs, ...
                            ind, val);
                    end
                otherwise
                    error('CHEBTECH:subsasgn:badType', ...
                        'Invalid subscripted assignment type.')
            end
        end

    end

end
