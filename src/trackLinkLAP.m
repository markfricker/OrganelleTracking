function tracks = trackLinkLAP(detections, p, wb, wbLo, wbHi)
%TRACKLINKLAP  Link per-frame detections into tracks using LAP assignment.
%
%   tracks = trackLinkLAP(detections, p)
%   tracks = trackLinkLAP(detections, p, wb, wbLo, wbHi)
%
%   Iterates over consecutive frame pairs, building and solving a LAP cost
%   matrix at each step.  Active track states (position, velocity) are
%   updated after each assignment.
%
%   Cost model (all terms normalised to [0,1] before weighting):
%     Displacement   - Euclidean distance to *predicted* position.
%                      Predicted = last position + smoothed velocity.
%                      Falls back to last position when the track is too
%                      short for reliable prediction.
%     Orientation    - Minimum angular difference of the major axis (mod pi).
%                      Handles the 180 deg ambiguity of regionprops orientation.
%     Shape          - |delta aspect-ratio| / p.maxAspectRatioNorm.
%     Area           - Relative area change (optional, toggle p.useAreaCost).
%     Direction      - Inconsistency between actual displacement and predicted
%                      velocity direction (optional, toggle p.useDirCost).
%
%   Spatial gate:
%     usePredictedGate = true  (default): gate is maxDisplacement from the
%       predicted position.  Tighter — appropriate for fast-moving objects.
%     usePredictedGate = false: legacy behaviour — gate is 2*maxDisplacement
%       from the last known position.  Wider, less strict.
%
%   Requires MATLAB Computer Vision Toolbox (assignDetectionsToTracks).
%
%   INPUTS
%     detections       : [nT x 1] cell array from trackDetectionsExtract.
%     p                : parameter struct (see trackParamsDefault).
%     wb, wbLo, wbHi  : optional waitbar handle and fractional range.
%                        Cancel signalled by setappdata(wb,'cancel',true).
%
%   OUTPUT
%     tracks : struct array, one element per track, with fields:
%       id            Unique integer track identifier.
%       frames        [1 x nObs] observed frame indices.
%       x, y          [1 x nObs] centroid positions (px).
%       vx, vy        [1 x nObs] velocity estimates (px/frame).
%                     NaN for the first observation of each track.
%       speed         [1 x nObs] ||v|| in px/frame (NaN first obs).
%       area          [1 x nObs]
%       majorAxis     [1 x nObs]
%       minorAxis     [1 x nObs]
%       aspectRatio   [1 x nObs]
%       orientation   [1 x nObs] degrees
%       eccentricity  [1 x nObs]
%       labelID       [1 x nObs] label in original stack

    useWb = nargin >= 3 && ~isempty(wb) && ishandle(wb);

    nT = numel(detections);

    % Internal track state: one row per active track.
    state  = struct('id',{},'lastX',{},'lastY',{},'vx',{},'vy',{}, ...
                    'nObs',{},'lastFrame',{});
    nextID = 1;

    % Collected output buffers -- cell per track, consolidated at end.
    buf = struct('id',{},'frames',{},'x',{},'y',{},'vx',{},'vy',{}, ...
                 'speed',{},'area',{},'majorAxis',{},'minorAxis',{}, ...
                 'aspectRatio',{},'orientation',{},'eccentricity',{}, ...
                 'labelID',{});

    maxDisp         = p.maxDisplacement;
    maxDispSq       = maxDisp^2;
    maxAR           = 10;           % normalisation ceiling for aspect ratio
    velAlpha        = p.velocityAlpha;
    velMinFrames    = p.velocityMinFrames;
    costNonAssign   = p.costNonAssignment;

    % Optional cost-term flags (default on if field missing).
    usePredGate = getf(p, 'usePredictedGate', true);
    useAreaCost = getf(p, 'useAreaCost',       true);
    useDirCost  = getf(p, 'useDirCost',        true);
    wArea       = getf(p, 'wArea',             0.2);
    wDir        = getf(p, 'wDir',              0.3);

    for iT = 1:nT

        if useWb
            frac = wbLo + (iT / nT) * (wbHi - wbLo);
            waitbar(frac, wb, sprintf('Linking frame %d / %d', iT, nT));
            if getappdata(wb, 'cancel')
                error('trackCompute:cancelled', ...
                    'Tracking cancelled by user at linking frame %d / %d.', iT, nT);
            end
        end

        det = detections{iT};

        %% Case: no detections this frame -- all active tracks persist.
        %  Gap handling is done later in trackGapClose; here we just skip.
        if isempty(det)
            continue
        end

        nDet    = numel(det.x);
        nTracks = numel(state);

        if nTracks == 0
            % Seed all detections as new tracks.
            for k = 1:nDet
                [state, buf, nextID] = createTrack(state, buf, nextID, ...
                    iT, det, k, NaN, NaN);
            end
            continue
        end

        %% Build cost matrix [nTracks x nDet].
        cost = inf(nTracks, nDet);

        % Pre-compute predicted positions for every active track.
        predX      = zeros(nTracks, 1);
        predY      = zeros(nTracks, 1);
        hasPred    = false(nTracks, 1);   % true when velocity estimate is reliable

        for i = 1:nTracks
            if state(i).nObs >= velMinFrames && ...
               ~(isnan(state(i).vx) || isnan(state(i).vy))
                predX(i)   = state(i).lastX + state(i).vx;
                predY(i)   = state(i).lastY + state(i).vy;
                hasPred(i) = true;
            else
                predX(i) = state(i).lastX;
                predY(i) = state(i).lastY;
            end
        end

        for i = 1:nTracks
            % Only link to tracks observed in the immediately preceding frame.
            % Tracks last seen further back are left for gap closing.
            if state(i).lastFrame < iT - 1
                continue
            end

            for j = 1:nDet

                % ---- Spatial gate -------------------------------------------
                if usePredGate
                    % Gate on predicted position with radius maxDisp.
                    dxP = det.x(j) - predX(i);
                    dyP = det.y(j) - predY(i);
                    if (dxP^2 + dyP^2) > maxDispSq
                        continue
                    end
                    cDisp = (dxP^2 + dyP^2) / maxDispSq;
                else
                    % Legacy: wider gate (2x) from last known position.
                    dx0 = det.x(j) - state(i).lastX;
                    dy0 = det.y(j) - state(i).lastY;
                    if (dx0^2 + dy0^2) > maxDispSq * 4
                        continue
                    end
                    % Displacement cost measured from predicted position.
                    dxP = det.x(j) - predX(i);
                    dyP = det.y(j) - predY(i);
                    cDisp = (dxP^2 + dyP^2) / maxDispSq;
                end

                % ---- Orientation cost (mod 180 deg) -------------------------
                dTheta = mod(abs(det.orientation(j) - ...
                              trackLastOrientation(buf, state(i).id)), 180);
                dTheta = min(dTheta, 180 - dTheta);   % range 0..90
                cOrient = (dTheta / 90)^2;

                % ---- Aspect-ratio cost --------------------------------------
                dAR    = abs(det.aspectRatio(j) - trackLastAR(buf, state(i).id));
                cShape = min(dAR / maxAR, 1)^2;

                % ---- Area cost (optional) -----------------------------------
                cArea = 0;
                if useAreaCost
                    lastA = trackLastArea(buf, state(i).id);
                    if lastA > 0
                        relDeltaA = abs(det.area(j) - lastA) / lastA;
                        cArea = min(relDeltaA, 1)^2;
                    end
                end

                % ---- Direction cost (optional) ------------------------------
                % Penalise when the candidate's displacement from the last
                % position points in a direction inconsistent with the track's
                % velocity.  Only applied when a reliable velocity exists.
                cDir = 0;
                if useDirCost && hasPred(i) && ...
                   ~(isnan(state(i).vx) || state(i).vx == 0 && state(i).vy == 0)
                    % Unit vector of track velocity.
                    vMag = hypot(state(i).vx, state(i).vy);
                    uvx  = state(i).vx / vMag;
                    uvy  = state(i).vy / vMag;
                    % Displacement from last known position.
                    dx0  = det.x(j) - state(i).lastX;
                    dy0  = det.y(j) - state(i).lastY;
                    dMag = hypot(dx0, dy0);
                    if dMag > 0
                        cosTheta = (dx0 * uvx + dy0 * uvy) / dMag;
                        % cosTheta in [-1,1]; 1 = perfect alignment.
                        % Map to cost: 0 when aligned, 1 when opposed.
                        cDir = ((1 - cosTheta) / 2)^2;
                    end
                end

                c = p.wDisp  * cDisp   + ...
                    p.wOrient * cOrient + ...
                    p.wShape  * cShape  + ...
                    wArea     * cArea   + ...
                    wDir      * cDir;

                % Final hard cutoff: normalised displacement > 1 is rejected.
                if cDisp > 1
                    c = Inf;
                end
                cost(i, j) = c;
            end
        end

        %% Solve LAP.
        % Guard: replace any NaN entries (e.g. from degenerate regionprops
        % properties on near-circular objects) with Inf so they are never
        % assigned.
        cost(isnan(cost)) = Inf;

        [assignments, unassignedTracks, unassignedDets] = ...
            assignDetectionsToTracks(cost, costNonAssign);

        %% Update assigned tracks.
        for r = 1:size(assignments, 1)
            ti = assignments(r, 1);
            di = assignments(r, 2);

            % Velocity update via EMA.
            newVx = det.x(di) - state(ti).lastX;
            newVy = det.y(di) - state(ti).lastY;
            if isnan(state(ti).vx)
                state(ti).vx = newVx;
                state(ti).vy = newVy;
            else
                state(ti).vx = velAlpha * state(ti).vx + (1-velAlpha) * newVx;
                state(ti).vy = velAlpha * state(ti).vy + (1-velAlpha) * newVy;
            end

            state(ti).lastX     = det.x(di);
            state(ti).lastY     = det.y(di);
            state(ti).nObs      = state(ti).nObs + 1;
            state(ti).lastFrame = iT;

            buf = appendDetection(buf, state(ti).id, iT, det, di, ...
                state(ti).vx, state(ti).vy);
        end

        %% Unassigned tracks: state persists unchanged (gap closing later).
        %  unassignedTracks intentionally unused here.  %#ok<NASGU>
        nUnassignedTracks = numel(unassignedTracks); %#ok<NASGU>

        %% Unassigned detections: start new tracks.
        for k = unassignedDets(:)'
            [state, buf, nextID] = createTrack(state, buf, nextID, ...
                iT, det, k, NaN, NaN);
        end
    end

    %% Consolidate buffers into output struct array.
    tracks = consolidateTracks(buf);
end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function val = trackLastOrientation(buf, id)
    idx = find([buf.id] == id, 1, 'last');
    if isempty(idx) || isempty(buf(idx).orientation)
        val = 0;
    else
        val = buf(idx).orientation(end);
    end
end

function val = trackLastAR(buf, id)
    idx = find([buf.id] == id, 1, 'last');
    if isempty(idx) || isempty(buf(idx).aspectRatio)
        val = 1;
    else
        val = buf(idx).aspectRatio(end);
    end
end

function val = trackLastArea(buf, id)
    idx = find([buf.id] == id, 1, 'last');
    if isempty(idx) || isempty(buf(idx).area)
        val = 0;
    else
        val = buf(idx).area(end);
    end
end

function [state, buf, nextID] = createTrack(state, buf, nextID, iT, det, k, vx, vy)
    ns.id        = nextID;
    ns.lastX     = det.x(k);
    ns.lastY     = det.y(k);
    ns.vx        = vx;
    ns.vy        = vy;
    ns.nObs      = 1;
    ns.lastFrame = iT;
    state(end+1) = ns; %#ok<AGROW>

    nb.id          = nextID;
    nb.frames      = iT;
    nb.x           = det.x(k);
    nb.y           = det.y(k);
    nb.vx          = NaN;
    nb.vy          = NaN;
    nb.speed       = NaN;
    nb.area        = det.area(k);
    nb.majorAxis   = det.majorAxis(k);
    nb.minorAxis   = det.minorAxis(k);
    nb.aspectRatio = det.aspectRatio(k);
    nb.orientation = det.orientation(k);
    nb.eccentricity= det.eccentricity(k);
    nb.labelID     = det.labelID(k);
    buf(end+1)     = nb; %#ok<AGROW>

    nextID = nextID + 1;
end

function buf = appendDetection(buf, id, iT, det, k, vx, vy)
    idx = find([buf.id] == id, 1, 'last');
    buf(idx).frames      = [buf(idx).frames,      iT];
    buf(idx).x           = [buf(idx).x,           det.x(k)];
    buf(idx).y           = [buf(idx).y,            det.y(k)];
    buf(idx).vx          = [buf(idx).vx,           vx];
    buf(idx).vy          = [buf(idx).vy,           vy];
    buf(idx).speed       = [buf(idx).speed,        hypot(vx, vy)];
    buf(idx).area        = [buf(idx).area,         det.area(k)];
    buf(idx).majorAxis   = [buf(idx).majorAxis,    det.majorAxis(k)];
    buf(idx).minorAxis   = [buf(idx).minorAxis,    det.minorAxis(k)];
    buf(idx).aspectRatio = [buf(idx).aspectRatio,  det.aspectRatio(k)];
    buf(idx).orientation = [buf(idx).orientation,  det.orientation(k)];
    buf(idx).eccentricity= [buf(idx).eccentricity, det.eccentricity(k)];
    buf(idx).labelID     = [buf(idx).labelID,      det.labelID(k)];
end

function tracks = consolidateTracks(buf)
    if isempty(buf)
        tracks = struct([]);
        return
    end
    tracks = buf;
    % Ensure row vectors.
    for i = 1:numel(tracks)
        fns = fieldnames(tracks(i));
        for f = 1:numel(fns)
            v = tracks(i).(fns{f});
            if isnumeric(v) && ~isscalar(v)
                tracks(i).(fns{f}) = v(:)';
            end
        end
    end
end

function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
