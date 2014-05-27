function a = truncate(a,tol)
% TRUNCATE(A), truncates off all elements in A which are below absolute
% value tol. 
    
            %fast binary search. 
            factor = 2;
            while(1)
                cut = floor(length(a)*(1-1/factor));
                if(cut<1 || length(a)<=1), break; end;
                if(max(abs(a(cut:end)))<tol)
                    a=a(1:cut);
                else
                    factor = factor*2;
                end
                if(factor>16), break;end
            end
            if(length(a)<=1), return; end
            % slower truncate to finish the job. 
            while(abs(a(end))<tol)
                a(end)=[];
            end

end