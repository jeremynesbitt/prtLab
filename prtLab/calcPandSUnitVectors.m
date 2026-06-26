function [p,s] = calcPandSUnitVectors(k, ada)
        % Deal with normal incidence
        if dot(k,ada) == 1 & abs(k - [0;0;1]) < eps 
            s = [1;0;0];
        else
        s = cross(k,ada);
        s = s/norm(s);
        end
        p = cross(k,s);
        p = p/norm(p);
end