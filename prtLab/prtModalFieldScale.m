function scale = prtModalFieldScale(modeE, fieldE)
%PRTMODALFIELDSCALE Least-squares scalar mapping a modal E field to fieldE.

denom = modeE' * modeE;
if abs(denom) < eps
    scale = 0;
else
    scale = (modeE' * fieldE) / denom;
end

end
