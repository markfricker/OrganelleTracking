function c = trackEvaluateLinkCost(baseTracks, labelStack, frame, srcLabel, tgtLabel, p)
%TRACKEVALUATELINKCOST  Score an arbitrary (frame,srcLabel)->(frame+1,tgtLabel)
%candidate using the SAME cost formula as trackLinkLAP, for diagnostic /
%learning purposes (e.g. comparing the cost of a manually-corrected link
%against the cost of whatever the tracker originally chose).
%
%   c = trackEvaluateLinkCost(baseTracks, labelStack, frame, srcLabel, tgtLabel, p)
%
%   baseTracks : struct array — the ORIGINAL (unedited) tracker output.
%                Must have per-observation vx/vy (the EMA-updated velocity
%                as computed live during the original trackLinkLAP run) —
%                reused here rather than re-deriving the full EMA chain.
%   labelStack : [nY nX nT] label matrix used for the original tracking run.
%   frame      : source frame (1-based). Target is looked up at frame+1.
%   srcLabel, tgtLabel : label values in labelStack.
%   p          : parameter struct (see trackParamsDefault) — must be the
%                SAME p used for the original run, so weights/flags match.
%
%   OUTPUT c : struct with fields
%     found        — true if srcLabel resolves to a track at `frame` and
%                    tgtLabel resolves to a real object at `frame+1`.
%     cDisp, cOrient, cShape, cArea, cDir, cBoundary — individual [0,1]
%                    cost terms (0 if the corresponding useXCost flag is off).
%     total        — weighted sum (Inf if cDisp > 1, matching trackLinkLAP's
%                    hard cutoff — i.e. this candidate would have been
%                    gated out entirely, never even reaching the LAP solve).
%     rejected     — logical, true iff total is Inf.
%
%   NOTE ON FIDELITY: this reconstructs the source track's pre-transition
%   state (position, velocity, last orientation/area/aspect ratio) from its
%   OWN stored observation at `frame` in baseTracks, rather than re-running
%   the sequential EMA/gating process. For a track observed continuously up
%   to `frame` this exactly reproduces the state trackLinkLAP had at
%   decision time. If `frame` falls inside a gap-closed (interpolated,
%   labelID==0) stretch, the reconstructed state instead reflects
%   trackGapClose's linear interpolation, which is a reasonable diagnostic
%   approximation but not a bit-exact replay of the original per-frame LAP.

c = struct('found', false, 'cDisp', NaN, 'cOrient', NaN, 'cShape', NaN, ...
    'cArea', NaN, 'cDir', NaN, 'cBoundary', NaN, 'total', NaN, 'rejected', false);

% ---- Resolve source track's state at `frame` ----------------------------
kSrc = 0;  idx = 0;
for k = 1:numel(baseTracks)
    fIdx = find(baseTracks(k).frames == frame, 1);
    if ~isempty(fIdx) && baseTracks(k).labelID(fIdx) == srcLabel
        kSrc = k; idx = fIdx;
        break
    end
end
if kSrc == 0, return; end
tr = baseTracks(kSrc);

lastX = tr.x(idx);  lastY = tr.y(idx);
vx    = tr.vx(idx); vy    = tr.vy(idx);
if isnan(vx), vx = 0; end
if isnan(vy), vy = 0; end
nObs  = idx;   % approximation: assumes no gap before `frame` (see docstring)

lastOrient = tr.orientation(idx);
lastAR     = tr.aspectRatio(idx);
lastArea   = tr.area(idx);

% ---- Resolve target detection's features at frame+1 ---------------------
if frame+1 > size(labelStack,3), return; end
mask = labelStack(:,:,frame+1) == tgtLabel;
if ~any(mask(:)), return; end
props = regionprops(mask, 'Centroid','Area','MajorAxisLength','MinorAxisLength','Orientation');
if isempty(props), return; end
pr  = props(1);
detX = pr.Centroid(1);  detY = pr.Centroid(2);
detArea = pr.Area;
detAR   = pr.MajorAxisLength / max(pr.MinorAxisLength, 1);
detOrient = pr.Orientation;
touchesBoundary = false;
if isfield(p,'boundaryLabels') && ~isempty(p.boundaryLabels) && ...
        frame+1 <= numel(p.boundaryLabels) && ~isempty(p.boundaryLabels{frame+1})
    touchesBoundary = ismember(tgtLabel, p.boundaryLabels{frame+1});
end

c.found = true;

% ---- Same weights/flags as trackLinkLAP ----------------------------------
maxDisp   = p.maxDisplacement;
maxDispSq = maxDisp^2;
maxAR     = 10;
velMinFrames = p.velocityMinFrames;
useAreaCost     = getf(p, 'useAreaCost',      true);
useDirCost      = getf(p, 'useDirCost',       true);
useBoundaryCost = getf(p, 'useBoundaryCost',  true);
wArea     = getf(p, 'wArea',      0.2);
wDir      = getf(p, 'wDir',       0.3);
wBoundary = getf(p, 'wBoundary',  0.3);

hasPred = nObs >= velMinFrames && ~(vx == 0 && vy == 0);
if hasPred
    predX = lastX + vx;
    predY = lastY + vy;
else
    predX = lastX;
    predY = lastY;
end

% ---- Displacement --------------------------------------------------------
dxP = detX - predX;  dyP = detY - predY;
c.cDisp = (dxP^2 + dyP^2) / maxDispSq;

% ---- Orientation ----------------------------------------------------------
dTheta = mod(abs(detOrient - lastOrient), 180);
dTheta = min(dTheta, 180 - dTheta);
c.cOrient = (dTheta / 90)^2;

% ---- Shape (aspect ratio) --------------------------------------------------
dAR = abs(detAR - lastAR);
c.cShape = min(dAR / maxAR, 1)^2;

% ---- Area (optional) -------------------------------------------------------
c.cArea = 0;
if useAreaCost && lastArea > 0
    relDeltaA = abs(detArea - lastArea) / lastArea;
    c.cArea = min(relDeltaA, 1)^2;
end

% ---- Direction (optional) --------------------------------------------------
c.cDir = 0;
if useDirCost && hasPred
    vMag = hypot(vx, vy);
    if vMag > 0
        uvx = vx / vMag;  uvy = vy / vMag;
        dx0 = detX - lastX;  dy0 = detY - lastY;
        dMag = hypot(dx0, dy0);
        if dMag > 0
            cosTheta = (dx0*uvx + dy0*uvy) / dMag;
            c.cDir = ((1 - cosTheta) / 2)^2;
        end
    end
end

% ---- Boundary (optional) ---------------------------------------------------
c.cBoundary = 0;
if useBoundaryCost && touchesBoundary
    c.cBoundary = 1;
end

% ---- Total ------------------------------------------------------------------
c.total = p.wDisp*c.cDisp + p.wOrient*c.cOrient + p.wShape*c.cShape + ...
    wArea*c.cArea + wDir*c.cDir + wBoundary*c.cBoundary;
if c.cDisp > 1
    c.total = Inf;
end
c.rejected = isinf(c.total);
end

function v = getf(s, field, default)
if isfield(s, field) && ~isempty(s.(field)), v = s.(field); else, v = default; end
end
