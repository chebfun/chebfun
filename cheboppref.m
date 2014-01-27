classdef cheboppref < chebpref
    
    % See above for documentation.
    properties
%         discretisation = @colloc2
%         scale = NaN
%         dimensionValues = [32 64 128 256 512 724 1024 1448 2048]
    end
    
    properties (Dependent)
        % TODO: Check how to get this working again.
        discretization
    end

    methods

        function outPref = cheboppref()           
            
            outPref = outPref@chebpref;
            
            % TODO: Check with AA about preferences
            outPref.prefList.maxTotalLength = 2500;
            outPref.prefList.enableSingularityDetection = false;  % not supported
            
            % Default new properties.
            outPref.prefList.discretisation = @colloc2;
            outPref.prefList.scale = NaN;
            outPref.prefList.dimensionValues =[32 64 128 256 512 724 1024 1448 2048];
            
            % Preferences for nonlinear ODEs
            outPref.prefList.damped = 1;
            outPref.prefList.display = 'iter';
            outPref.prefList.errTol = 1e-10;
            outPref.prefList.lambdaMin = 1e-6;
            outPref.prefList.maxIter = 25;
            outPref.prefList.plotting = 1;
        end
        
        function out = get.discretization(pref)
            out = pref.discretisation;
        end

        function pref = set.discretization(pref,disc)
            pref.discretisation = disc;
        end

        function out = subsref(pref, ind)
        %SUBSREF   Subscripted referencing for CHEBOPPREF.
        %   P.PROP, where P is a CHEBOPPREF object, returns the value of the
        %   CHEBOPPREF property PROP stored in P.  If PROP is not a CHEBOPPREF
        %   property, P.TECHPREFS.PROP will be returned instead.  If PROP is
        %   neither a CHEBOPPREF property nor a field in P.TECHPREFS, an error
        %   will be thrown.
        %
        %   For access to fields PROP of TECHPREFS that have the same name as a
        %   CHEBOPPREF property, use the syntax P.TECHPREFS.PROP.
        %
        %   CHEBOPPREF does not support any other subscripted referencing types,
        %   including '()' and '{}'.
            switch ( ind(1).type )
                case '.'
                    if ( isfield(pref.prefList, ind(1).subs) )
                        out = pref.prefList.(ind(1).subs);
                    else
                        out = pref.prefList.techPrefs.(ind(1).subs);
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
        %SUBSASGN   Subscripted assignment for CHEBOPPREF.
        %   P.PROP = VAL, where P is a CHEBOPPREF object, assigns the value VAL
        %   to the CHEBOPPREF property PROP stored in P.  If PROP is not a
        %   CHEBOPPREF property, the assignment will be made to P.TECHPREFS.PROP
        %   instead.
        %
        %   To assign to fields PROP of TECHPREFS that have the same name as a
        %   CHEBOPPREF property, use the syntax P.TECHPREFS.PROP = VAL.
        %
        %   CHEBOPPREF does not support any other subscripted assignment types,
        %   including '()' and '{}'.
            switch ( ind(1).type )
                case '.'
                    if ( isfield(pref.prefList, ind(1).subs) )
                        pref.prefList = builtin('subsasgn', pref.prefList, ...
                            ind, val);
                    else
                        pref.prefList.techPrefs = builtin('subsasgn', ...
                            pref.prefList.techPrefs, ind, val);
                    end
                otherwise
                    error('CHEBTECH:subsasgn:badType', ...
                        'Invalid subscripted assignment type.')
            end
        end

    end

end
