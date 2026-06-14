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
%   Invariant tracks (spatially identical in ALL five configurations) are
%   drawn first as faint grey lines so they provide spatial context without
%   dominating.  Variable tracks are drawn on top in colour, with 'x'
%   markers at speed-spike frames.  This focuses attention on the objects
%   where the toggles actually have an effect.
%
%   Colours are assigned from the "all on" reference so the same physical
%   object has the same colour across columns.  Tracks with no match in the
%   reference (only present in this config) are shown in bright red.
%
%   Title of each column: config name + "nTr tracks / nVar variable / nSpk spikes"
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

    OVERLAP_THR  = 0.70;
    DIST_THR     = 3.0;
    INV_ALPHA    = 0.18;   % opacity of invariant (background) lines
    INV_COLOUR   = [0.55 0.55 0.55];
    NOVEL_COLOUR = [0.9 0.15 0.15];   % tracks present here but not in reference

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

    % ---- Classify each track as invariant or variable ----------------------
    isVar = cell(1, nCfg);
    for k = 1:nCfg
        nTr      = numel(allTracks{k});
        isVar{k} = true(1, nTr);
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
                isVar{k}(t) = false;
            end
        end
    end

    % ---- Colour map from reference ------------------------------------------
    refTracks = allTracks{REF};
    nRef      = numel(refTracks);
    cmap      = lines(max(nRef, 1));
    if nRef > size(cmap,1)
        cmap = repmat(cmap, ceil(nRef/size(cmap,1)), 1);
        cmap = cmap(1:nRef,:);
    end

    % ---- Kymograph background -----------------------------------------------
    if ~isempty(imStack) && isequal(size(imStack), [nY nX nT])
        kymo = squeeze(max(double(imStack), [], 1))';   % [nT x nX]
        kymo = mat2gray(kymo);
    else
        kymo = 0.15 * ones(nT, nX);
    end

    % ---- Plot — single row of kymographs ------------------------------------
    fig = figure('Name','Toggle Comparison','NumberTitle','off', ...
                 'Position',[30 30 300*nCfg 480]);
    tl  = tiledlayout(fig, 1, nCfg, 'TileSpacing','compact','Padding','compact');
    title(tl, 'LAP tracker — cost toggle comparison', 'FontWeight','bold');

    axK = gobjects(1, nCfg);

    for k = 1:nCfg
        tracks = allTracks{k};
        nTr    = numel(tracks);
        invIdx = find(~isVar{k});
        varIdx = find( isVar{k});
        nVar   = numel(varIdx);
        nSpk   = countSpikes(tracks(varIdx));

        axK(k) = nexttile(tl, k);
        imagesc(axK(k), kymo);
        colormap(axK(k), gray);
        axis(axK(k), 'tight');
        hold(axK(k), 'on');

        % ---- Invariant tracks first (faint grey, no spikes) -----------------
        for ti = invIdx
            tr = tracks(ti);
            plot(axK(k), tr.x, tr.frames, '-', ...
                 'Color', [INV_COLOUR, INV_ALPHA], 'LineWidth', 1);
        end

        % ---- Variable tracks on top (coloured, spikes marked) ---------------
        for ti = varIdx
            tr  = tracks(ti);
            col = matchColour(tr, refTracks, cmap, NOVEL_COLOUR, OVERLAP_THR, DIST_THR);
            plot(axK(k), tr.x, tr.frames, '-', 'Color', col, 'LineWidth', 1.5);

            spikeF = spikeFrames(tr);
            if ~isempty(spikeF)
                [~, idx] = ismember(spikeF, tr.frames);
                idx = idx(idx > 0);
                plot(axK(k), tr.x(idx), tr.frames(idx), 'x', ...
                     'Color', col, 'MarkerSize', 9, 'LineWidth', 2);
            end
        end

        ttl = sprintf('%s\n%d tracks / %d var / %d spk', ...
                      configs{k,4}, nTr, nVar, nSpk);
        title(axK(k), ttl, 'FontWeight','bold', 'FontSize', 9);
        if k == 1
            ylabel(axK(k), 'Frame');
        else
            set(axK(k), 'YTickLabel', {});
        end
        xlabel(axK(k), 'X (px)');
    end

    linkaxes(axK, 'xy');
    fprintf('Done.\n');
end

% -------------------------------------------------------------------------
% Matching helpers
% -------------------------------------------------------------------------

function tf = anyMatch(tr, candidates, overlapThr, distThr)
    tf = false;
    for k = 1:numel(candidates)
        if tracksMatch(tr, candidates(k), overlapThr, distThr)
            tf = true;
            return
        end
    end
end

function tf = tracksMatch(a, b, overlapThr, distThr)
    sharedF = intersect(a.frames, b.frames);
    nShared = numel(sharedF);
    if nShared == 0, tf = false; return, end
    if nShared / min(numel(a.frames), numel(b.frames)) < overlapThr
        tf = false; return
    end
    [~, ia] = ismember(sharedF, a.frames);
    [~, ib] = ismember(sharedF, b.frames);
    tf = mean(hypot(a.x(ia) - b.x(ib), a.y(ia) - b.y(ib))) < distThr;
end

function col = matchColour(tr, refTracks, cmap, novelCol, overlapThr, distThr)
    for k = 1:numel(refTracks)
        if tracksMatch(tr, refTracks(k), overlapThr, distThr)
            col = cmap(mod(k-1, size(cmap,1))+1, :);
            return
        end
    end
    col = novelCol;
end

function frames = spikeFrames(tr)
    spd = tr.speed;
    spd(isnan(spd)) = 0;
    frames = [];
    if numel(spd) < 3, return, end
    med = median(spd(spd > 0));
    if med == 0, return, end
    spikeObs = find(abs(diff(spd)) > 2*med) + 1;
    frames   = tr.frames(spikeObs);
end

function n = countSpikes(tracks)
    n = 0;
    for k = 1:numel(tracks)
        if ~isempty(spikeFrames(tracks(k))), n = n + 1; end
    end
end
