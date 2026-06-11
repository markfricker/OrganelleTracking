function mitoSimPlayMovie(frames, groundTruth, showGroundTruth)
%MITOSIMPLAYMOVIE  Display a synthetic mitochondria movie with optional overlay.
%
%   mitoSimPlayMovie(frames)
%   mitoSimPlayMovie(frames, groundTruth)
%   mitoSimPlayMovie(frames, groundTruth, showGroundTruth)
%
%   INPUTS
%     frames          : [nY nX nT] uint16 image stack (from mitoSimGenerate).
%     groundTruth     : [nT x 1] cell array of mito struct arrays (optional).
%                       Pass [] to omit overlay.
%     showGroundTruth : logical -- whether to draw the GT overlay (default true).
%
%   Layout: two panels side-by-side.  Left = raw image, right = GT overlay.
%   Red ellipses show the true mito boundary; red cross marks centroid; ID
%   label printed in red.

    if nargin < 2, groundTruth     = []; end
    if nargin < 3, showGroundTruth = true; end

    nT       = size(frames, 3);
    hasGT    = showGroundTruth && ~isempty(groundTruth);
    nPanels  = 1 + hasGT;

    hFig = figure('Name', 'Synthetic Mitochondria Movie', ...
                  'Position', [50, 50, 500*nPanels, 520]);

    for t = 1:nT
        if ~ishandle(hFig), break, end

        subplot(1, nPanels, 1);
        imshow(frames(:,:,t), []);
        title(sprintf('Frame %d / %d', t, nT));

        if hasGT
            subplot(1, nPanels, 2);
            imshow(frames(:,:,t), []);
            hold on;

            mito = groundTruth{t};
            for i = 1:numel(mito)
                drawEllipse(mito(i));
            end
            hold off;
            title('Ground Truth Overlay');
        end

        drawnow;
        pause(0.05);
    end
end

% -------------------------------------------------------------------------

function drawEllipse(m)
    % Draw a rotated ellipse outline + centroid marker + ID label.
    phi = linspace(0, 2*pi, 60);
    xe  = m.majorAxis .* cos(phi);
    ye  = m.minorAxis .* sin(phi);

    xr  =  xe .* cos(m.theta) - ye .* sin(m.theta) + m.x;
    yr  =  xe .* sin(m.theta) + ye .* cos(m.theta) + m.y;

    plot(xr, yr, 'r-', 'LineWidth', 1);
    plot(m.x, m.y, 'r+', 'MarkerSize', 8, 'LineWidth', 1.5);
    text(m.x + 5, m.y - 5, num2str(m.id), ...
        'Color', 'red', 'FontSize', 8);
end
