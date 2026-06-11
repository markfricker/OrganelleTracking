function demoTrackMethods(labelStack, imRaw, p)
%DEMOTRACKMETHODS  Visualise tracking results and diagnose streaming events.
%
%   demoTrackMethods(labelStack, imRaw, p)
%
%   Produces a 3-panel figure:
%     Panel 1  Trajectory overlay on the mean-intensity projection of imRaw.
%              Each track is drawn as a coloured trail; streaming segments
%              (speed > p.streamingThreshold) are highlighted in red.
%              Track ID labels at the final position.
%     Panel 2  Speed (µm/s) vs frame for every track; streaming segments
%              are marked with red circles.  Grey band = mean ± 1 SD of all
%              instantaneous speeds (diffusive regime reference).
%     Panel 3  Orientation (°) vs frame; streaming tracks highlighted.
%
%   INPUTS
%     labelStack : [nY nX nT] uint16 from segmentation.
%     imRaw      : [nY nX nT] numeric raw intensity (any class).
%                  Pass labelStack if raw image unavailable — trajectory
%                  overlay will be shown on the label image instead.
%     p          : parameter struct (see trackParamsDefault).
%                  Call with p = trackParamsDefault() as a starting point.
%
%   EXAMPLE
%     p = trackParamsDefault();
%     p.minTrackLength = 3;
%     demoTrackMethods(myLabels, myRaw, p);

    if nargin < 3 || isempty(p)
        p = trackParamsDefault();
    end
    if nargin < 2 || isempty(imRaw)
        imRaw = labelStack;
    end

    %% Run tracking.
    t0 = tic;
    tracks = trackCompute(labelStack, p);
    elapsed = toc(t0);

    nT = size(labelStack, 3);
    fprintf('trackCompute: %d tracks from %d frames in %.1f s\n', ...
        numel(tracks), nT, elapsed);

    if isempty(tracks)
        warning('demoTrackMethods: no tracks produced.');
        return
    end

    %% Background image: mean projection, normalised.
    bg = mat2gray(mean(double(imRaw), 3));

    %% Colour map: one colour per track (cycling through lines colormap).
    nTracks = numel(tracks);
    cmap    = lines(max(nTracks, 1));
    cmap    = cmap(mod((1:nTracks)-1, size(cmap,1))+1, :);

    %% ---- Panel 1: trajectory overlay ----
    fig = figure('Name', 'Tracking — trajectory & speed overview', ...
                 'NumberTitle', 'off', ...
                 'Color', 'k', ...
                 'Units', 'normalized', ...
                 'OuterPosition', [0.02 0.02 0.96 0.96]);

    t = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(t);
    imshow(bg, [], 'Parent', ax1);
    hold(ax1, 'on');
    axis(ax1, 'image', 'off');
    title(ax1, sprintf('Trajectories  n=%d tracks', nTracks), ...
          'Color', 'w', 'FontSize', 9);

    for i = 1:nTracks
        tr  = tracks(i);
        col = cmap(i, :);

        % Draw full trail as thin coloured line.
        plot(ax1, tr.x, tr.y, '-', 'Color', [col, 0.5], 'LineWidth', 0.8);

        % Highlight streaming segments.
        if any(tr.streaming)
            xs = tr.x;  xs(~tr.streaming) = NaN;
            ys = tr.y;  ys(~tr.streaming) = NaN;
            plot(ax1, xs, ys, 'r-', 'LineWidth', 2);
        end

        % Label at final position.
        text(ax1, tr.x(end), tr.y(end), sprintf('%d', tr.id), ...
             'Color', col, 'FontSize', 6, 'HorizontalAlignment', 'center');
    end
    hold(ax1, 'off');

    %% ---- Panel 2: speed vs frame ----
    ax2 = nexttile(t);
    hold(ax2, 'on');
    set(ax2, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    title(ax2, 'Speed (µm/s)', 'Color', 'w', 'FontSize', 9);
    xlabel(ax2, 'Frame', 'Color', 'w');
    ylabel(ax2, 'Speed (µm/s)', 'Color', 'w');

    % Collect all instantaneous speeds for reference band.
    allSpeeds = [];
    for i = 1:nTracks
        sp = tracks(i).speedUm;
        allSpeeds = [allSpeeds, sp(~isnan(sp))]; %#ok<AGROW>
    end
    if ~isempty(allSpeeds)
        mu  = mean(allSpeeds);
        sg  = std(allSpeeds);
        xlim_val = [1, nT];
        patch(ax2, [xlim_val(1), xlim_val(2), xlim_val(2), xlim_val(1)], ...
              [mu-sg, mu-sg, mu+sg, mu+sg], ...
              [0.3 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        yline(ax2, mu, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
    end

    for i = 1:nTracks
        tr  = tracks(i);
        col = cmap(i, :);
        plot(ax2, tr.frames, tr.speedUm, '-', 'Color', [col, 0.6], 'LineWidth', 0.8);
        % Mark streaming frames.
        if any(tr.streaming)
            plot(ax2, tr.frames(tr.streaming), tr.speedUm(tr.streaming), ...
                 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
        end
    end
    hold(ax2, 'off');

    %% ---- Panel 3: orientation vs frame ----
    ax3 = nexttile(t);
    hold(ax3, 'on');
    set(ax3, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    title(ax3, 'Orientation (°)', 'Color', 'w', 'FontSize', 9);
    xlabel(ax3, 'Frame', 'Color', 'w');
    ylabel(ax3, 'Orientation (°)', 'Color', 'w');
    ylim(ax3, [-90, 90]);

    for i = 1:nTracks
        tr  = tracks(i);
        col = cmap(i, :);
        plot(ax3, tr.frames, tr.orientation, '-', 'Color', [col, 0.6], 'LineWidth', 0.8);
        if any(tr.streaming)
            plot(ax3, tr.frames(tr.streaming), tr.orientation(tr.streaming), ...
                 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
        end
    end
    hold(ax3, 'off');

    linkaxes([ax2, ax3], 'x');
    xlim(ax2, [1, nT]);
end
