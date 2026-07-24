function interaction = traceSurfaceInteraction(T, surfaceIndex, ray, hit, normal, options)
%TRACESURFACEINTERACTION Dispatch one ray/surface interaction.

mediumIn = T(surfaceIndex,:);
mediumOut = T(surfaceIndex+1,:);
caseName = classifyInterfaceCase(mediumIn.MaterialType, mediumOut.MaterialType);

interaction = emptyInteractionRecord();
interaction.surfaceIndex = surfaceIndex;
interaction.caseName = caseName;
interaction.position = hit;
interaction.normal = normal;
interaction.incident = makeIncidentRecord(ray);
interaction.frames = makeFrameRecord();

switch caseName
    case "isotropicToIsotropic"
        interaction = traceIsotropicToIsotropic( ...
            interaction, mediumIn, mediumOut, ray, hit, normal, options);

    case "isotropicToUniaxial"
        interaction = traceIsotropicToUniaxial( ...
            interaction, mediumIn, mediumOut, ray, hit, normal, options);

    case "uniaxialToIsotropic"
        interaction = traceUniaxialToIsotropic( ...
            interaction, mediumIn, mediumOut, ray, hit, normal, options);

    case "uniaxialToUniaxial"
        interaction = traceUniaxialToUniaxial( ...
            interaction, mediumIn, mediumOut, ray, hit, normal, options);

    otherwise
        error('traceSurfaceInteraction:UnsupportedCase', ...
            'Interface case "%s" is not supported.', caseName);
end

interaction.exiting = makeExitingRecords(interaction.children);

function incident = makeIncidentRecord(ray)
incident = struct();
incident.rayId = ray.id;
incident.mode = ray.mode;
incident.branchType = ray.branchType;
incident.mediumType = ray.mediumType;
incident.position = ray.position;
incident.k = ray.k;
incident.S = ray.S;
incident.modeE = ray.modeE;
incident.modeH = ray.modeH;
incident.fieldE = ray.fieldE;
incident.fieldH = ray.fieldH;
incident.E = ray.E;
incident.H = ray.H;
incident.P = ray.P;
incident.Q = ray.Q;
incident.O = ray.O;
incident.localBasis = ray.localBasis;
incident.amplitude = ray.amplitude;
incident.flux = ray.flux;
incident.OPL = ray.OPL;
end

function frames = makeFrameRecord()
frames = struct();
frames.Oin = [];
frames.Oout = [];
frames.Qin = [];
frames.Qout = [];
frames.inputBasis = struct();
frames.outputBasis = struct();
end

function exiting = makeExitingRecords(children)
exiting = repmat(struct( ...
    'mode', "", ...
    'branchType', "", ...
    'mediumType', "", ...
    'position', zeros(3,1), ...
    'k', zeros(3,1), ...
    'S', zeros(3,1), ...
    'modeE', zeros(3,1), ...
    'modeH', zeros(3,1), ...
    'fieldE', zeros(3,1), ...
    'fieldH', zeros(3,1), ...
    'E', zeros(3,1), ...
    'H', zeros(3,1), ...
    'P', eye(3), ...
    'Q', eye(3), ...
    'O', eye(3), ...
    'localBasis', struct(), ...
    'amplitude', [], ...
    'flux', [], ...
    'OPL', []), 0, 1);

for jj = 1:numel(children)
    exiting(jj).mode = children(jj).mode;
    exiting(jj).branchType = children(jj).branchType;
    exiting(jj).mediumType = children(jj).mediumType;
    exiting(jj).position = children(jj).position;
    exiting(jj).k = children(jj).k;
    exiting(jj).S = children(jj).S;
    exiting(jj).modeE = children(jj).modeE;
    exiting(jj).modeH = children(jj).modeH;
    exiting(jj).fieldE = children(jj).fieldE;
    exiting(jj).fieldH = children(jj).fieldH;
    exiting(jj).E = children(jj).E;
    exiting(jj).H = children(jj).H;
    exiting(jj).P = children(jj).P;
    exiting(jj).Q = children(jj).Q;
    exiting(jj).O = children(jj).O;
    exiting(jj).localBasis = children(jj).localBasis;
    exiting(jj).amplitude = children(jj).amplitude;
    exiting(jj).flux = children(jj).flux;
    exiting(jj).OPL = children(jj).OPL;
end
end

end
