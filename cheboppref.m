classdef cheboppref < chebpref
    
    % See above for documentation.
    properties
        damped
        display
        errTol
        lambdaMin
        maxIter
        plotting
        discretisation = @colloc2
        scale = NaN
        dimensionValues = [32 64 128 256 512 724 1024 1448 2048]
    end
    
    properties (Dependent)
        discretization
    end

    methods

        function outPref = cheboppref()           
            
            outPref = outPref@chebpref;
            
            outPref.maxTotalLength = 2500;
            outPref.enableSingularityDetection = false;  % not supported
        end
        
        function out = get.discretization(pref)
            out = pref.discretisation;
        end

        function pref = set.discretization(pref,disc)
            pref.discretisation = disc;
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
