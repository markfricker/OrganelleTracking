function trackMergeMontage(labelStack, mergeInfo, tracks, imRaw, p)
%TRACKMERGEMONTAGE  Show each resolved transient merge as a before|split|after strip.
%
%   trackMergeMontage(labelStack, mergeInfo, tracks, imRaw)
%   trackMergeMontage(labelStack, mergeInfo, tracks, imRaw, p)
%
%   For every split recorded in mergeInfo (from trackResolveMerges /
%   trackCompute with p.resolveMerges), builds ONE concatenated image per
%   merge — the frame before, the resolved frame, and the frame after,
%   cropped to the same window and placed side by side.  Object outlines are
%   drawn over the raw image; at the split frame the host object is outlined
%   CYAN and the newly split-off object MAGENTA, others grey.  The host and
%   new objects' centroids are joined by a trajectory line ACROSS the three
%   panels (cyan / magenta), so you can see the two objects converge into the
%   merge and separate again.
%
%   INPUTS
%     labelStack : [nY nX nT] uint16 — the merge-resolved stack.
%     mergeInfo  : struct from trackCompute with field .splits.
%     tracks     : struct array from trackCompute (for centroid trajectories).
%     imRaw      : [nY nX nT] backdrop image, or [] to use the mask.
%     p          : optional; p.montageMax caps rows (default 10).
%
%   See also: trackResolveMerges, trackCompute, trackOutlineProjection

    if nargin < 5, p = struct(); end
    montageMax = getf(p, 'montageMax', 24);

    if ~isfield(mergeInfo, 'splits') || isempty(mergeInfo.splits)
        fprintf('trackMergeMontage: no resolved merges to show.\n');
        return
    end

    [nY, nX, nT] = size(labelStack);
    splits = mergeInfo.splits;
    nShow  = min(numel(splits), montageMax);

    % Grid layout: each tile is one merge's before|split|after strip.
    nCols = max(1, round(getf(p, 'montageCols', ceil(sqrt(nShow)))));
    nRows = ceil(nShow / nCols);

    colHost = [0 1 1];      % cyan   (host centroid trajectory)
    colNew  = [1 0 1];      % magenta (new  centroid trajectory)
    colEdge = [1 1 0];      % yellow  (object boundaries)

    figTitle = sprintf('Resolved merges  (%d of %d shown)', nShow, numel(splits));
    figure('Name', figTitle, 'NumberTitle','off', 'Color','k', ...
           'Units','normalized', 'OuterPosition',[0.04 0.05 0.92 0.9]);
    tl = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
    title(tl, figTitle, 'Color','w', 'FontWeight','bold', 'Interpreter','none');

    for s = 1:nShow
        sp = splits(s);
        f  = sp.frame;

        % Shared crop window around the split pair (+ pad).
        m = (labelStack(:,:,f) == sp.hostLabel) | (labelStack(:,:,f) == sp.newLabel);
        [rr, cc] = find(m);
        if isempty(rr), continue; end
        pad = 15;
        r0 = max(1, min(rr)-pad);  r1 = min(nY, max(rr)+pad);
        c0 = max(1, min(cc)-pad);  c1 = min(nX, max(cc)+pad);
        H = r1-r0+1;  W = c1-c0+1;

        gFrames = [f-1, f, f+1];

        % Greyscale backdrop strip (panels concatenated side by side).
        strip = zeros(H, 3*W, 'single');
        for k = 1:3
            g = gFrames(k);
            if g < 1 || g > nT, continue; end
            off = (k-1)*W;
            if ~isempty(imRaw)
                bg = mat2gray(single(imRaw(r0:r1, c0:c1, g)));
            else
                bg = single(labelStack(r0:r1, c0:c1, g) > 0) * 0.5;
            end
            strip(:, off+(1:W)) = bg;
        end

        ax = nexttile(tl);
        imshow(strip, [], 'Parent', ax); hold(ax, 'on');

        % Object boundaries as closed yellow line plots (per label, so the
        % internal split line is drawn, not just the outer union).
        for k = 1:3
            g = gFrames(k);
            if g < 1 || g > nT, continue; end
            off = (k-1)*W;
            Lc  = labelStack(r0:r1, c0:c1, g);
            for lab = unique(Lc(Lc > 0))'
                B = bwboundaries(Lc == lab, 'noholes');
                for b = 1:numel(B)
                    bd = B{b};
                    plot(ax, bd(:,2)+off, bd(:,1), '-', 'Color', colEdge, 'LineWidth', 1.5);
                end
            end
        end

        % panel separators + labels
        for k = 1:2, xline(ax, k*W + 0.5, 'Color',[0.4 0.4 0.4], 'LineWidth',0.5); end
        if s <= nCols
            lbl = {'before','split','after'};
            for k = 1:3
                text(ax, (k-1)*W + W/2, 3, lbl{k}, 'Color','w', 'FontSize',9, ...
                     'HorizontalAlignment','center', 'VerticalAlignment','top');
            end
        end
        ylabel(ax, sprintf('f%d: L%d\\rightarrowL%d+L%d', f, sp.hostLabel, sp.hostLabel, sp.newLabel), ...
               'Color','y', 'FontSize',9);
        set(ax, 'YColor','k', 'XColor','k', 'Visible','on', 'XTick',[], 'YTick',[]);

        % centroid trajectories across the three panels
        drawTrajectory(ax, findTrackAt(tracks, f, sp.hostLabel), gFrames, r0, c0, W, colHost);
        drawTrajectory(ax, findTrackAt(tracks, f, sp.newLabel),  gFrames, r0, c0, W, colNew);

        hold(ax, 'off');
    end
end

% -------------------------------------------------------------------------

function drawTrajectory(ax, tr, gFrames, r0, c0, W, col)
%DRAWTRAJECTORY  Join a track's centroid across the three time panels.
    if isempty(tr), return; end
    pts = [];
    for k = 1:3
        j = find(tr.frames == gFrames(k), 1);
        if isempty(j), continue; end
        x = tr.x(j) - c0 + 1 + (k-1)*W;
        y = tr.y(j) - r0 + 1;
        pts = [pts; x y]; %#ok<AGROW>
    end
    if isempty(pts), return; end
    plot(ax, pts(:,1), pts(:,2), '-o', 'Color', col, ...
         'MarkerSize', 4, 'LineWidth', 1.2, 'MarkerFaceColor', col);
end

function tr = findTrackAt(tracks, frame, label)
%FINDTRACKAT  Track passing through (frame, label), or [] if none.
    tr = [];
    for t = 1:numel(tracks)
        j = find(tracks(t).frames == frame, 1);
        if ~isempty(j) && tracks(t).labelID(j) == label
            tr = tracks(t);
            return
        end
    end
end

function v = getf(s, field, default)
    if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
