function trackOutlineProjection(labelStack, tracks, imRaw, p)
%TRACKOUTLINEPROJECTION  Whole-movie projection of object outlines by track.
%
%   trackOutlineProjection(labelStack, tracks, imRaw)
%   trackOutlineProjection(labelStack, tracks, imRaw, p)
%
%   Flattens the entire movie into a single image: every tracked object's
%   OUTLINE at every frame, coloured by its track ID, painted over a
%   max-projection backdrop, with the centroid trajectory of each track
%   drawn as a line in the same colour.  Gives a one-glance overview of
%   where every object went — including objects recovered by merge
%   resolution (their split frames carry the new track's colour).
%
%   INPUTS
%     labelStack : [nY nX nT] uint16 — the (merge-resolved) label stack.
%     tracks     : struct array from trackCompute (fields frames, x, y,
%                  labelID; labelID==0 = gap-interpolated, no outline drawn).
%     imRaw      : [nY nX nT] backdrop image, or [] for a black backdrop.
%     p          : optional struct:
%                    outlineEveryN  paint outline every Nth frame (default 1).
%                    fadeWithTime   dim earlier frames so motion direction
%                                   reads as dark->bright (default true).
%                    showTrails     overlay centroid trajectory lines (default true).
%
%   See also: trackMergeMontage, trackMovie, trackVerify

    if nargin < 4, p = struct(); end
    everyN = max(1, round(getf(p, 'outlineEveryN', 1)));
    fade   = getf(p, 'fadeWithTime', true);
    trails = getf(p, 'showTrails',   true);

    [nY, nX, nT] = size(labelStack);
    nTr = numel(tracks);
    if nTr == 0
        fprintf('trackOutlineProjection: no tracks to show.\n');
        return
    end

    % Backdrop: dim max-projection of the raw image (or black).
    if ~isempty(imRaw)
        bg = mat2gray(max(single(imRaw), [], 3)) * 0.4;
    else
        bg = zeros(nY, nX, 'single');
    end
    R = bg; G = bg; B = bg;

    % Distinct colour per track (shuffled hsv so neighbours differ).
    C = hsv(max(nTr, 1));
    rng(0); C = C(randperm(nTr), :);

    % Paint outlines, max-blended so overlaps keep the brighter colour.
    for t = 1:nTr
        tr  = tracks(t);
        col = C(t, :);
        for i = 1:numel(tr.frames)
            L = tr.labelID(i);
            f = tr.frames(i);
            if L == 0 || mod(f-1, everyN) ~= 0
                continue
            end
            per = bwperim(labelStack(:, :, f) == L);
            if ~any(per(:)), continue; end
            w = 1;
            if fade, w = 0.35 + 0.65 * (f / nT); end
            R(per) = max(R(per), col(1) * w);
            G(per) = max(G(per), col(2) * w);
            B(per) = max(B(per), col(3) * w);
        end
    end
    RGB = cat(3, R, G, B);

    figTitle = sprintf('Outline projection — %d tracks, %d frames', nTr, nT);
    figure('Name', figTitle, 'NumberTitle','off', 'Color','k', ...
           'Units','normalized', 'OuterPosition',[0.1 0.1 0.7 0.8]);
    ax = axes('Parent', gcf);
    imshow(RGB, 'Parent', ax); hold(ax, 'on');
    title(ax, figTitle, 'Color','w', 'FontWeight','bold', 'Interpreter','none');

    if trails
        for t = 1:nTr
            tr = tracks(t);
            plot(ax, tr.x, tr.y, '-', 'Color', [C(t,:) 0.5], 'LineWidth', 1);
            % start (o) and end (s) markers
            plot(ax, tr.x(1),   tr.y(1),   'o', 'Color', C(t,:), 'MarkerSize', 4, 'LineWidth', 1);
            plot(ax, tr.x(end), tr.y(end), 's', 'Color', C(t,:), 'MarkerSize', 5, 'LineWidth', 1);
        end
    end
    hold(ax, 'off');
end

% -------------------------------------------------------------------------
function v = getf(s, field, default)
    if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
