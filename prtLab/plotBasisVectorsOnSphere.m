function plotBasisVectorsOnSphere(k_grid, x_loc, y_loc, opts)
% plotBasisVectorsOnSphere - Plot local basis vectors on unit propagation sphere
%
% plotBasisVectorsOnSphere(k_grid, x_loc, y_loc) plots x (red) and y (green)
% basis vectors as arrows at each point k on the unit sphere.
%
% plotBasisVectorsOnSphere(k_grid, x_loc, [], opts) plots only x-basis vectors.
%
% Inputs:
%   k_grid - Nx3 array of unit propagation vectors (points on sphere)
%   x_loc  - Nx3 array of local x-basis vectors
%   y_loc  - Nx3 array of local y-basis vectors (pass [] to skip)
%   opts   - struct with options:
%     .showSphere  (default true)   - draw transparent unit sphere
%     .xColor      (default 'r')    - color for x-basis vectors
%     .yColor      (default [0 0.7 0]) - color for y-basis vectors
%     .arrowScale  (default 0.15)   - length of arrows
%     .sphereAlpha (default 0.15)   - sphere face transparency
%     .sphereColor (default [0.7 0.7 1]) - sphere face color
%     .sphereN     (default 30)     - sphere mesh resolution
%     .viewAngle   (default [-37.5, 30]) - view angle [az, el]
%     .lineWidth   (default 1)      - arrow line width

if nargin < 4, opts = struct(); end
if ~isfield(opts, 'showSphere'),  opts.showSphere = true; end
if ~isfield(opts, 'xColor'),      opts.xColor = 'r'; end
if ~isfield(opts, 'yColor'),      opts.yColor = [0 0.7 0]; end
if ~isfield(opts, 'arrowScale'),  opts.arrowScale = 0.15; end
if ~isfield(opts, 'sphereAlpha'), opts.sphereAlpha = 0.15; end
if ~isfield(opts, 'sphereColor'), opts.sphereColor = [0.7 0.7 1]; end
if ~isfield(opts, 'sphereN'),     opts.sphereN = 30; end
if ~isfield(opts, 'viewAngle'),   opts.viewAngle = [-37.5, 30]; end
if ~isfield(opts, 'lineWidth'),   opts.lineWidth = 1; end

% Draw unit sphere
if opts.showSphere
    [X, Y, Z] = sphere(opts.sphereN);
    surf(X, Y, Z, 'FaceAlpha', opts.sphereAlpha, 'EdgeAlpha', 0.08, ...
        'FaceColor', opts.sphereColor);
    hold on;
end

s = opts.arrowScale;

% Plot x-basis vectors
if ~isempty(x_loc)
    quiver3(k_grid(:,1), k_grid(:,2), k_grid(:,3), ...
        x_loc(:,1)*s, x_loc(:,2)*s, x_loc(:,3)*s, 0, ...
        'Color', opts.xColor, 'LineWidth', opts.lineWidth, 'MaxHeadSize', 0.3);
    hold on;
end

% Plot y-basis vectors
if ~isempty(y_loc)
    quiver3(k_grid(:,1), k_grid(:,2), k_grid(:,3), ...
        y_loc(:,1)*s, y_loc(:,2)*s, y_loc(:,3)*s, 0, ...
        'Color', opts.yColor, 'LineWidth', opts.lineWidth, 'MaxHeadSize', 0.3);
end

xlabel('x'); ylabel('y'); zlabel('z');
axis equal;
view(opts.viewAngle);
grid on;

end
