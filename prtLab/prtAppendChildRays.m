function result = prtAppendChildRays(rays, parentRayId, interaction, options, childMetadata)
%PRTAPPENDCHILDRAYS Filter and register active interaction children.

arguments
    rays struct
    parentRayId (1,1) double
    interaction (1,1) struct
    options (1,1) struct
    childMetadata (1,1) struct = struct()
end

childIds = [];
metadataNames = fieldnames(childMetadata);

for ii = 1:numel(interaction.children)
    child = interaction.children(ii);
    if ~child.active || child.flux < options.minFlux || ...
            abs(child.amplitude) < options.minAmplitude
        continue;
    end

    child.id = numel(rays) + 1;
    child.parentId = parentRayId;
    child.history = [rays(parentRayId).history, child.id];
    for jj = 1:numel(metadataNames)
        name = metadataNames{jj};
        child.metadata.(name) = childMetadata.(name);
    end
    rays(end+1) = child; %#ok<AGROW>
    childIds(end+1) = child.id; %#ok<AGROW>
end

result = struct('rays', rays, 'ids', childIds);

end
