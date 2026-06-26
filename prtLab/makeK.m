function K = makeK(k_inc)

K = [0, -k_inc(3), k_inc(2); ...
        k_inc(3), 0, -k_inc(1) ; ...
        -k_inc(2), k_inc(1), 0];

end