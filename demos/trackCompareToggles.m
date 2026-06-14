function trackCompareToggles(labelStack, imStack, p)
%TRACKCOMPARETOGGLES  Compare LAP tracking with each cost toggle on/off.
%
%   trackCompareToggles(labelStack, imStack, p)
%
%   Re-runs trackCompute five times on the same labelStack, varying only
%   the three toggle flags (usePredictedGate, useAreaCost, useDirCost):
%
%     Col 1  All off       — baseline, legacy 2x-gate behaviour
%     Col 2  +PredGate     — tighter gate from predicted position only
%     Col 3  +AreaCost     — adds area-change penalty
%     Col 4  +DirCost      — adds motion-direction penalty
%     Col 5  All on        — full model (current defaults)
%
%   Tracks that are spatially identical across ALL five configurations are
%   suppressed — only tracks that change under at least one toggle are shown.
%   This focuses the comparison on the objects where the toggles actually
%   have an effect.
%
%   Each column shows a kymograph (max-proj along Y) with variable track
%   centroids overlaid, coloured consistently by their ID in the "all on"
%   reference.  A summary line below each column shows:
%     nTr    = total tracks (surviving minTrackLength filter)
%     nVar   = number of variable tracks shown
%     medL   = median track length (frames)
%     spikes = tracks with ≥1 speed jump > 2× their own median speed
%
%   INPUTS
%     labelStack : [nY nX nT] uint16 label matrix from organelle segment.
%     imStack    : [nY nX nT] raw image for kymograph background.
%                  Pass [] to use a flat grey background.
%     p          : parameter struct (see trackParamsDefault).  Toggle
%                  fields are overridden per column; all other fields kept.
%
%   MATCHING PARAMETERS (edit here to adjust sensitivity)
%     OVERLAP_THR  Minimum fraction of frames in common for a match [0,1].
%     DIST_THR     Maximum mean centroid distance on shared frames (px).

    OVERLAP_THR = 0.70;
    DIST_THR    = 3.0;

    if nargin < 2, imStack = []; end
    if nargin < 3, p = trackParamsDefault(); end

    configs = { ...
        false, false, false, 'All off';    ...
        true,  false, false, '+PredGate';  ...
        false, true,  false, '+AreaCost';  ...
        false, false, true,  '+DirCost';   ...
        true,  true,  true,  'All on';     ...
    };
    nCfg = size(configs, 1);
    REF  = nCfg;   % "All on" is the reference for colour assignment

    [nY, nX, nT] = size(labelStack);

    % ---- Run tracker for each configuration --------------------------------
    fprintf('Running %d configurations...\n', nCfg);
    allTracks = cell(1, nCfg);
    for k = 1:nCfg
        pCfg                  = p;
        pCfg.usePredictedGate = configs{k,1};
        pCfg.useAreaCost      = configs{k,2};
        pCfg.useDirCost       = configs{k,3};
        allTracks{k}          = trackCompute(labelStack, pCfg);
        fprintf('  %s: %d tracks\n', configs{k,4}, numel(allTracks{k}));
    end

    % ---- Find invariant tracks using the reference (all on) ----------------
    % A track in config k is invariant if it matches a track in EVERY other config.
    % We build a per-config logical mask: isVar(k){t} = true if track t is variable.
    isVar = cell(1, nCfg);
    for k = 1:nCfg
        nTr       = numel(allTracks{k});
        isVar{k}  = true(1, nTr);   % start as all variable
        for t = 1:nTr
            matchedAll = true;
            for kk = 1:nCfg
                if kk == k, continue, end
                if ~anyMatch(allTracks{k}(t), allTracks{kk}, OVERLAP_THR, DIST_THR)
                    matchedAll = false;
                    break
                end
            end
            if matchedAll
                isVar{k}(t) = false;   % identical in all configs — suppress
            end
        end
    end

    % ---- Assign colours from reference config --------------------------------
    % Tracks in other configs are coloured by their best match in the reference.
    % Unmatched tracks (new in this config) get a fixed grey.
    refTracks = allTracks{REF};
    nRef      = numel(refTracks);
    cmap      = lines(max(nRef, 1));
    if nRef > size(cmap,1)
        cmap = repmat(cmap, ceil(nRef/size(cmap,1)), 1);
        cmap = cmap(1:nRef,:);
    end
    greyCol = [0.6 0.6 0.6];

    % ---- Build kymograph background -----------------------------------------
    if ~isempty(imStack) && isequal(size(imStack), [nY nX nT])
        kymo = squeeze(max(double(imStack), [], 1))';   % [nT x nX]
        kymo = mat2gray(kymo);
    else
        kymo = 0.15 * ones(nT, nX);
    end

    % ---- Plot ---------------------------------------------------------------
    fig = figure('Name','Toggle Comparison','NumberTitle','off', ...
                 'Position',[30 30 300*nCfg 520]);
    tl  = tiledlayout(fig, 2, nCfg, 'TileSpacing','compact','Padding','compact');
    title(tl, 'LAP tracker — cost toggle comparison (variable tracks only)', ...
          'FontWeight','bold');

    axK = gobjects(1, nCfg);

    for k = 1:nCfg
        tracks = allTracks{k};
        nTr    = numel(tracks);
        varIdx = find(isVar{k});
        nVar   = numel(varIdx);

        % Summary stats (all tracks, not just variable).
        if nTr > 0
            lengths = arrayfun(@(t) numel(t.frames), tracks);
            medL    = median(lengths);
            spikes  = countSpikes(tracks);
        else
            medL = 0;  spikes = 0;
        end

        % ---- Kymograph tile (top row) ----------------------------------------
        axK(k) = nexttile(tl, k);
        imagesc(axK(k), kymo);
        colormap(axK(k), gray);
        axis(axK(k), 'tight');
        hold(axK(k), 'on');

        for ti = varIdx
            tr   = tracks(ti);
            col  = matchColour(tr, refTracks, cmap, greyCol, OVERLAP_THR, DIST_THR);

            % Line joining centroids along time axis.
            plot(axK(k), tr.x, tr.frames, '-', 'Color', col, 'LineWidth', 1.5);

            % Spike markers: frames where speed jumps > 2× local median.
            spikeF = spikeFrames(tr);
            if ~isempty(spikeF)
                [~, idx] = ismember(spikeF, tr.frames);
                idx = idx(idx > 0);
                plot(axK(k), tr.x(idx), tr.frames(idx), 'x', ...
                     'Color', col, 'MarkerSize', 9, 'LineWidth', 2);
            end
        end

        title(axK(k), configs{k,4}, 'FontWeight','bold');
        if k == 1
            ylabel(axK(k), 'Frame');
        else
            set(axK(k), 'YTickLabel', {});
        end
        xlabel(axK(k), 'X (px)');

        % ---- Stats tile (bottom row) -----------------------------------------
        nexttile(tl, nCfg + k);
        axis off
        txt = sprintf('tracks: %d\nshown:  %d\nmedLen: %.0f fr\nspikes: %d', ...
                      nTr, nVar, medL, spikes);
        text(0.5, 0.5, txt, 'Units','normalized', ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'FontSize', 10);
    end

    linkaxes(axK, 'xy');
    fprintf('Done. Invariant tracks suppressed; showing only variable tracks.\n');
end

% -------------------------------------------------------------------------
% Matching helpers
% -------------------------------------------------------------------------

function tf = anyMatch(tr, candidates, overlapThr, distThr)
%ANYMATCH  True if tr matches at least one track in candidates.
    tf = false;
    for k = 1:numel(candidates)
        if tracksMatch(tr, candidates(k), overlapThr, distThr)
            tf = true;
            return
        end
    end
end

function tf = tracksMatch(a, b, overlapThr, distThr)
%TRACKSMATCH  True if tracks a and b share enough frames with close centroids.
    sharedF = intersect(a.frames, b.frames);
    nShared = numel(sharedF);
    if nShared == 0
        tf = false;
        return
    end
    % Fraction of the shorter track covered by shared frames.
    overlapFrac = nShared / min(numel(a.frames), numel(b.frames));
    if overlapFrac < overlapThr
        tf = false;
        return
    end
    % Mean centroid distance on shared frames.
    [~, ia] = ismember(sharedF, a.frames);
    [~, ib] = ismember(sharedF, b.frames);
    dx   = a.x(ia) - b.x(ib);
    dy   = a.y(ia) - b.y(ib);
    dist = mean(hypot(dx, dy));
    tf   = dist < distThr;
end

function col = matchColour(tr, refTracks, cmap, greyCol, overlapThr, distThr)
%MATCHCOLOUR  Return colour from refTracks colour map, or grey if unmatched.
    for k = 1:numel(refTracks)
        if tracksMatch(tr, refTracks(k), overlapThr, distThr)
            col = cmap(mod(k-1, size(cmap,1))+1, :);
            return
        end
    end
    col = greyCol;
end

% -------------------------------------------------------------------------
function frames = spikeFrames(tr)
%SPIKEFRAMES  Return frame numbers where speed jumps by > 2× local median.
%   The jump at observation i is assigned to frame i+1 (the later frame).
    spd = tr.speed;
    spd(isnan(spd)) = 0;
    frames = [];
    if numel(spd) < 3, return, end
    med = median(spd(spd > 0));
    if med == 0, return, end
    jumps = abs(diff(spd));          % length nObs-1; jump(i) is between obs i and i+1
    spikeObs = find(jumps > 2*med) + 1;   % index into tr.frames
    frames = tr.frames(spikeObs);
end

function n = countSpikes(tracks)
%COUNTSPIKES  Count tracks with ≥1 speed spike.
    n = 0;
    for k = 1:numel(tracks)
        if ~isempty(spikeFrames(tracks(k)))
            n = n + 1;
        end
    end
end
