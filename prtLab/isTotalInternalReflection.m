function tf = isTotalInternalReflection(k, normal, n1, n2)
cosThetaI = dot(k, normal);
sin2ThetaI = max(0, 1 - cosThetaI^2);
tf = isreal(n1) && isreal(n2) && n1 > n2 && (n1/n2)^2 * sin2ThetaI > 1;
end
