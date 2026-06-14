function interaction = traceUniaxialToUniaxialStub(interaction, mediumIn, mediumOut, ray, hit, normal, options)
%TRACEUNIAXIALTOUNIAXIALSTUB Placeholder for uniaxial-to-uniaxial tracing.

if options.dispatchUnimplemented == "error"
    error('traceUniaxialToUniaxialStub:NotImplemented', ...
        'uniaxial-to-uniaxial tracing is not implemented yet.');
end

childO = makeChildTemplate(ray, hit, normal, mediumOut, "ordinary");
childE = makeChildTemplate(ray, hit, normal, mediumOut, "extraordinary");
interaction.children = [childO; childE];

end
