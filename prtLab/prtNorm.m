function q = prtNorm(v)
%prtNorm normalization of vectors, whether real or complex
q = v/sqrt(sum(v.*v));
end