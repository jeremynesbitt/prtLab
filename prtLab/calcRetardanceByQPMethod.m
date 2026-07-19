function retardance = calcRetardanceByQPMethod(Q,P, method)
% From PLAOS Chapter 17
% method; hopefully temp var to eval how retardance is calculated
% optional parameter, 1 will subtract before calculating angle

[U,~,V] = svd(Q\P);

eigVals = eig(V\U);

x = 10*eps;

% One of these eigenvalues is 1.  Need to find it
realVals = real(eigVals);
fazVals  = angle(eigVals);
if abs(abs(realVals(1))-1) < x 
    retVals = eigVals(2:3);
elseif abs(abs(realVals(2))-1) < x 
    retVals = [eigVals(1);eigVals(3)];
elseif abs(abs(realVals(3))-1) < x 
    retVals = eigVals(1:2);
end

if exist('retVals', 'var') == 0 
    disp('Warning!  Could not find Eigenvalue=1 within eps.  Loosen requirement')
    x = 1E-5;
    if abs(abs(realVals(1))-1) < x 
        retVals = eigVals(2:3);
    elseif abs(abs(realVals(2))-1) < x 
        retVals = [eigVals(1);eigVals(3)];
    elseif abs(abs(realVals(3))-1) < x 
        retVals = eigVals(1:2);
    end   
end



% if abs(abs(realVals(1))-1) < x && abs(fazVals(1)) < x
%     retVals = eigVals(2:3);
% elseif abs(abs(realVals(2))-1) < x && abs(fazVals(2)) < x
%     retVals = [eigVals(1);eigVals(3)];
% elseif abs(abs(realVals(3))-1) < x && abs(fazVals(3)) < x
%     retVals = eigVals(1:2);
% end

% Enforce postive retardance
%eigVals = eig(inv(Q)*P);
%eigVals = eig(Q\P);

a1 = angle(retVals(1));
a2 = angle(retVals(2));

if a2 < a1
    retardance = a1-a2;
    if(exist('method', 'var'))
        if method==1
            retardance = angle(retVals(1)-retVals(2));
        end
    end

else
    retardance = a2-a1;
    if(exist('method', 'var'))
        if method==1
            retardance = angle(retVals(2)-retVals(1));
        end
    end    
end


%retardance = angle(eigVals(2))-angle(eigVals(3));

end