function area = domainvolume(f)
%DOMAINVOLUME    Volume of the domain of f
%
%   DOMAINVOLUME(F) returns the volume of the topological domain of f.

if ( isempty(f) )
    area = 0;
else
    dom = f.domain; 
    area = diff(dom(1:2)) * diff(dom(3:4)) * diff(dom(5:6));
end

end