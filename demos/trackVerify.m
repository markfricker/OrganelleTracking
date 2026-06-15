function trackVerify(tracks, imStack, p)
%TRACKVERIFY  Diagnostic plots for verifying LAP tracking correctness.
%
%   trackVerify(tracks)
%   trackVerify(tracks, imStack)
%   trackVerify(tracks, imStack, p)
%
%   Produces four panels designed to make identity swaps, fragmentation,
%   and erroneous long jumps visually obvious without manually scrubbing
%   every frame of the movie:
%
%   Panel 1 — Position-vs-frame spaghetti (X and Y separately)
%     Each track is one coloured line.  Identity swaps appear as sharp
%     discontinuities or crossed trajectories; premature termination
%     shows up as short line segments in a region where long tracks exist.
%
%   Panel 2 — Speed-vs-frame profile per track
%     Plots ||v|| (px/frame) vs frame for every track, coloured by track
%     ID.  Spurious links typically produce isolated speed spikes; correct
%     streaming tracks have smooth, sustained elevation.
%
%   Panel 3 — Kymograph (space-time slice)
%     A maximum-intensity projection along Y, so each streaming mitochondrion
%     appears as a diagonal streak.  Track centroids are overlaid as coloured
%     dots.  Slope = speed; parallelism = same velocity class; an abrupt
%     colour switch across a streak indicates an identity swap.
%     Requires imStack ([nY nX nT] double or uint16).
%
%   Panel 4 — Displacement jump histogram
%     Distribution of frame-to-frame displacements across all tracks.
%     A heavy tail or a secondary mode beyond the expected maxDisplacement
%     suggests wrong links or missed detections.
%
%   INPUTS
%     tracks   : struct array from trackLinkLAP / trackGapClose / trackFilter.
%     imStack  : optional [nY nX nT] image stack.  Panel 3 skipped if absent.
%     p        : optional parameter struct.  Uses p.maxDisplacement to draw
%                a reference line on Panel 4.
%
%   TIPS FOR READING THE PLOTS
%     - X/Y spaghetti: look for sharp kinks or tracks that suddenly change
%       gradient.  Continuous curves with similar gradient = correct stream.
%     - Speed profile: jumps to high speed in an otherwise slow track, or
%       single-frame speed drops in a fast track, are swap signatures.
%     - Kymograph: streaks should be single-coloured along their length.
%       Colour transitions mid-streak mean the track switched identity.
%     - Displacement histogram: the bulk should sit well below maxDisp.
%       A secondary peak near maxDisp means the gate is too generous or
%       objects are genuinely fast — increase maxDisplacement accordingly.

    if nargin < 2, imStack = []; end
    if nargin < 3, p       = struct(); end

    maxDisp = getf(p, 'maxDisplacement', []);

    nT = 0;
    for k = 1:numel(tracks)
        nT = max(nT, max(tracks(k).frames));
    end

    % Colour palette — one colour per track, recycled.
    nTr  = numel(tracks);
    cmap = lines(max(nTr, 1));
    if nTr > size(cmap,1)
        cmap = repmat(cmap, ceil(nTr/size(cmap,1)), 1);
        cmap = cmap(1:nTr, :);
    end

    hasStack = ~isempty(imStack) && ndims(imStack) == 3; %#ok<ISMAT>

    fig = figure('Name','Track Verification','NumberTitle','off', ...
                 'Position',[50 50 1400 900]);
    tl  = tiledlayout(fig, 2, 3, 'TileSpacing','compact','Padding','compact');

    % -------------------------------------------------------------------------
    % Panel 1: X-position and Y-position vs frame
    % -------------------------------------------------------------------------
    ax1 = nexttile(tl);
    hold(ax1,'on');
    ax2 = nexttile(tl);
    hold(ax2,'on');

    for k = 1:nTr
        tr  = tracks(k);
        col = cmap(mod(k-1, size(cmap,1))+1, :);
        conf = getConfirmed(tr);
        plotTrackLine(ax1, tr.frames, tr.x, col, conf, 1);
        plotTrackLine(ax2, tr.frames, tr.y, col, conf, 1);
    end
    xlabel(ax1,'Frame');  ylabel(ax1,'X (px)');  title(ax1,'X-position vs frame');
    xlabel(ax2,'Frame');  ylabel(ax2,'Y (px)');  title(ax2,'Y-position vs frame');
    xlim(ax1,[1 nT]);  xlim(ax2,[1 nT]);
    grid(ax1,'on');     grid(ax2,'on');

    % -------------------------------------------------------------------------
    % Panel 2: Speed profile per track
    % -------------------------------------------------------------------------
    ax3 = nexttile(tl);
    hold(ax3,'on');

    allSpeeds = [];
    for k = 1:nTr
        tr  = tracks(k);
        spd = tr.speed;
        spd(isnan(spd)) = 0;
        allSpeeds = [allSpeeds, spd]; %#ok<AGROW>
        col  = cmap(mod(k-1, size(cmap,1))+1, :);
        conf = getConfirmed(tr);
        plotTrackLine(ax3, tr.frames, spd, col, conf, 1.5);

        % Mark spike frames with a cross.
        spikeF = spikeFrames(tr);
        if ~isempty(spikeF)
            [~, idx] = ismember(spikeF, tr.frames);
            idx = idx(idx > 0);
            plot(ax3, tr.frames(idx), spd(idx), 'o', ...
                 'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 12, 'LineWidth', 1);
            plot(ax3, tr.frames(idx), spd(idx), 'x', ...
                 'Color', col, 'MarkerSize', 9, 'LineWidth', 2);
        end
    end
    xlabel(ax3,'Frame');  ylabel(ax3,'Speed (px/frame)');
    title(ax3,'Per-track speed profile');
    if ~isempty(maxDisp)
        yline(ax3, maxDisp, 'r--', 'maxDisp');
    end
    xlim(ax3,[1 nT]);
    grid(ax3,'on');

    % -------------------------------------------------------------------------
    % Panel 3: Displacement-jump histogram
    % -------------------------------------------------------------------------
    ax4 = nexttile(tl);

    jumps = [];
    for k = 1:nTr
        tr = tracks(k);
        dx = diff(tr.x);
        dy = diff(tr.y);
        jumps = [jumps, hypot(dx, dy)]; %#ok<AGROW>
    end
    if ~isempty(jumps)
        histogram(ax4, jumps, 30, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none');
        if ~isempty(maxDisp)
            xline(ax4, maxDisp, 'r--', sprintf('maxDisp=%g', maxDisp));
        end
    end
    xlabel(ax4,'Frame-to-frame displacement (px)');
    ylabel(ax4,'Count');
    title(ax4,'Displacement jump histogram');
    grid(ax4,'on');

    % -------------------------------------------------------------------------
    % Panel 4: Kymograph (max-proj along Y, centroids overlaid)
    % -------------------------------------------------------------------------
    if hasStack
        ax5 = nexttile(tl);

        % Max-intensity projection along Y → [nX x nT] image.
        kymo = squeeze(max(double(imStack), [], 1))';   % [nT x nX]
        imagesc(ax5, kymo);
        colormap(ax5, gray);
        axis(ax5,'tight');
        hold(ax5,'on');

        for k = 1:nTr
            tr  = tracks(k);
            col = cmap(mod(k-1, size(cmap,1))+1, :);
            plotTrackLine(ax5, tr.x, tr.frames, col, ...
                          getConfirmed(tr), 1.5);
            spikeF = spikeFrames(tr);
            if ~isempty(spikeF)
                [~, idx] = ismember(spikeF, tr.frames);
                idx = idx(idx > 0);
                plot(ax5, tr.x(idx), tr.frames(idx), 'o', ...
                     'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 12, 'LineWidth', 1);
                plot(ax5, tr.x(idx), tr.frames(idx), 'x', ...
                     'Color', col, 'MarkerSize', 9, 'LineWidth', 2);
            end
        end
        xlabel(ax5,'X (px)');  ylabel(ax5,'Frame');
        title(ax5,'Kymograph (max-proj Y) + track centroids');
    end

    title(tl, sprintf('Track Verification  —  %d tracks, %d frames', nTr, nT));
end

function conf = getConfirmed(tr)
%GETCONFIRMED  Return confirmed vector, or all-true if field absent.
    if isfield(tr, 'confirmed')
        conf = tr.confirmed;
    else
        conf = true(1, numel(tr.frames));
    end
end

function plotTrackLine(ax, xData, yData, col, conf, lw)
%PLOTTRACKLINE  Plot a track as a continuous solid line with diamond markers
%   at unconfirmed link junctions (fwd-bwd disagrees at that frame).
%   This keeps the track visually whole while flagging uncertain joints.
    plot(ax, xData, yData, '-', 'Color', col, 'LineWidth', lw);
    if ~all(conf)
        unconf = find(~conf);
        plot(ax, xData(unconf), yData(unconf), 'd', ...
             'Color', col, 'MarkerSize', 6, 'LineWidth', 1.2, ...
             'MarkerFaceColor', 'w');
    end
end

function frames = spikeFrames(tr)
    spd = tr.speed;
    spd(isnan(spd)) = 0;
    frames = [];
    if numel(spd) < 3, return, end
    med = median(spd(spd > 0));
    if med == 0, return, end
    jumps    = abs(diff(spd));
    spikeObs = find(jumps > 2*med) + 1;
    frames   = tr.frames(spikeObs);
end

function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
