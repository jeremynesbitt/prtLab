function k_t = snell_vector(k_i, eta, n1, n2)
% Vector form of Snell's law (Eq. 10.20 from Chipman)
% k_i: incident unit vector, eta: surface normal into transmitted medium (unit)
% n1: incident index, n2: transmitted index
    cos_i = dot(k_i, eta);
    nr = n1/n2;
    k_t = nr*k_i - (nr*cos_i - sqrt(1 - nr^2*(1 - cos_i^2)))*eta;
    k_t = k_t / norm(k_t);
end
