function tracks = trackGapClose(tracks, detections, p)
%TRACKGAPCLOSE  Second-pass LAP gap closing over track segment endpoints.
%
%   tracks = trackGapClose(tracks, detections, p)
%
%   After frame-to-frame linking, track segments can be fragmented at
%   frames where a mitochondrion is temporarily invisible (bleaching, brief
%   out-of-focus movement, occlusion).  This function does a global second
%   pass: it pairs each track-segment *end* with compatible segment *starts*
%   within p.gapMax frames, solves a second LAP, and merges matched pairs.
%
%   Cost model
%     Displacement   - Euclidean distance between end position and the
%                      position linearly extrapolated from the track's
%                      terminal velocity, evaluated at the start frame of the
%                      candidate.  Scaled by gap length x p.gapCostScale.
%     Orientation    - Same angular metric as trackLinkLAP.
%
%   Any gap that would require displacement > p.maxDisplacement x gap x 1.5
%   is forbidden (Inf cost).
%
%   INPUTS
%     tracks     : struct array from trackLinkLAP.
%     detections : [nT x 1] cell from trackDetectionsExtract (reserved for
%                  future per-frame orientation look-up at boundaries).
%     p          : parameter struct (see trackParamsDefault).
%
%   OUTPUT
%     tracks : struct array with merged segments and short tracks filtered
%              by p.minTrackLength.
%
%   Note: merged tracks inherit the *earlier* segment's id.  The later
%   segment is removed from the array.  Observations from the gap frames
%   are interpolated linearly and flagged with labelID = 0.

    if isempty(tracks)
        return
    end

    % detections reserved for future per-frame orientation look-up
    nTracks = numel(tracks);
    gapMax  = p.gapMax;
    maxDisp = p.maxDisplacement;
    scale   = p.gapCostScale;

    % Build index: endFrame and startFrame for each track segment.
    endFrames   = arrayfun(@(t) t.frames(end),  tracks);
    startFrames = arrayfun(@(t) t.frames(1),    tracks);

    % Identify candidate (end, start) pairs within gapMax.
    % Rows = ends, cols = starts.
    nEnds   = nTracks;
    nStarts = nTracks;
    cost    = inf(nEnds, nStarts);

    for e = 1:nEnds
        for s = 1:nStarts
            if s == e, continue, end
            gap = startFrames(s) - endFrames(e);
            if gap < 1 || gap > gapMax + 1
                continue
            end

            % Extrapolate end track position to start frame.
            vx = tracks(e).vx(end);
            vy = tracks(e).vy(end);
            if isnan(vx), vx = 0; end
            if isnan(vy), vy = 0; end
            predX = tracks(e).x(end) + vx * gap;
            predY = tracks(e).y(end) + vy * gap;

            dx = tracks(s).x(1) - predX;
            dy = tracks(s).y(1) - predY;
            dist = hypot(dx, dy);

            maxAllowed = maxDisp * gap * 1.5;
            if dist > maxAllowed
                continue
            end

            cDisp   = (dist / maxDisp)^2;
            dTheta  = mod(abs(tracks(s).orientation(1) - tracks(e).orientation(end)), 180);
            dTheta  = min(dTheta, 180 - dTheta);
            cOrient = (dTheta / 90)^2;

            cost(e, s) = (p.wDisp * cDisp + p.wOrient * cOrient) * (scale ^ (gap - 1));
        end
    end

    % Solve LAP -- treat each track end/start as a "detection"/"track".
    % Guard against NaN (e.g. from degenerate regionprops orientation values).
    cost(isnan(cost)) = Inf;
    if ~all(isinf(cost(:)))
        [assignments, ~, ~] = assignDetectionsToTracks(cost, p.costNonAssignment);
    else
        assignments = zeros(0, 2);
    end

    % Merge matched pairs (earlier segment absorbs later segment).
    toDelete = false(nTracks, 1);
    for r = 1:size(assignments, 1)
        e = assignments(r, 1);   % end track (earlier)
        s = assignments(r, 2);   % start track (later)
        if toDelete(e) || toDelete(s), continue, end

        gap = startFrames(s) - endFrames(e);

        % Interpolate gap frames.
        [gapFrames, gapX, gapY, gapVx, gapVy] = interpolateGap( ...
            tracks(e), tracks(s), gap);

        % Append gap + start-track observations to end track.
        tracks(e) = appendSegment(tracks(e), gapFrames, gapX, gapY, ...
            gapVx, gapVy, tracks(s));

        toDelete(s) = true;
    end

    tracks(toDelete) = [];

    % Filter by minimum track length.
    lengths = arrayfun(@(t) numel(t.frames), tracks);
    tracks(lengths < p.minTrackLength) = [];
end

% -------------------------------------------------------------------------

function [gf, gx, gy, gvx, gvy] = interpolateGap(tEnd, tStart, gap)
    % Linear interpolation of position across the gap.
    nGap = gap - 1;
    if nGap == 0
        gf = []; gx = []; gy = []; gvx = []; gvy = [];
        return
    end
    frac  = (1:nGap) / gap;
    gf    = tEnd.frames(end)  + (1:nGap);
    gx    = tEnd.x(end)   + frac * (tStart.x(1)   - tEnd.x(end));
    gy    = tEnd.y(end)   + frac * (tStart.y(1)    - tEnd.y(end));
    gvx   = repmat(tEnd.vx(end), 1, nGap);
    gvy   = repmat(tEnd.vy(end), 1, nGap);
end

function tEnd = appendSegment(tEnd, gf, gx, gy, gvx, gvy, tStart)
    % Gap frames (labelID = 0 = interpolated).
    nGap = numel(gf);
    if nGap > 0
        tEnd.frames      = [tEnd.frames,      gf];
        tEnd.x           = [tEnd.x,           gx];
        tEnd.y           = [tEnd.y,            gy];
        tEnd.vx          = [tEnd.vx,           gvx];
        tEnd.vy          = [tEnd.vy,           gvy];
        tEnd.speed       = [tEnd.speed,        hypot(gvx, gvy)];
        tEnd.area        = [tEnd.area,         repmat(mean(tEnd.area),      1, nGap)];
        tEnd.majorAxis   = [tEnd.majorAxis,    repmat(mean(tEnd.majorAxis), 1, nGap)];
        tEnd.minorAxis   = [tEnd.minorAxis,    repmat(mean(tEnd.minorAxis), 1, nGap)];
        tEnd.aspectRatio = [tEnd.aspectRatio,  repmat(mean(tEnd.aspectRatio),1,nGap)];
        tEnd.orientation = [tEnd.orientation,  repmat(tEnd.orientation(end),1, nGap)];
        tEnd.eccentricity= [tEnd.eccentricity, repmat(mean(tEnd.eccentricity),1,nGap)];
        tEnd.labelID     = [tEnd.labelID,      zeros(1, nGap, 'uint16')];
    end
    % Append start-segment observations.
    tEnd.frames      = [tEnd.frames,      tStart.frames];
    tEnd.x           = [tEnd.x,           tStart.x];
    tEnd.y           = [tEnd.y,            tStart.y];
    tEnd.vx          = [tEnd.vx,           tStart.vx];
    tEnd.vy          = [tEnd.vy,           tStart.vy];
    tEnd.speed       = [tEnd.speed,        tStart.speed];
    tEnd.area        = [tEnd.area,         tStart.area];
    tEnd.majorAxis   = [tEnd.majorAxis,    tStart.majorAxis];
    tEnd.minorAxis   = [tEnd.minorAxis,    tStart.minorAxis];
    tEnd.aspectRatio = [tEnd.aspectRatio,  tStart.aspectRatio];
    tEnd.orientation = [tEnd.orientation,  tStart.orientation];
    tEnd.eccentricity= [tEnd.eccentricity, tStart.eccentricity];
    tEnd.labelID     = [tEnd.labelID,      tStart.labelID];
end
