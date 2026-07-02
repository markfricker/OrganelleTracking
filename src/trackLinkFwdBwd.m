function [tracks, confirmedKeys] = trackLinkFwdBwd(detections, p, wb, wbLo, wbHi, pins)
%TRACKLINKFWDBWD  Forward-backward LAP linking with confirmation scoring.
%
%   [tracks, confirmedKeys] = trackLinkFwdBwd(detections, p)
%   [tracks, confirmedKeys] = trackLinkFwdBwd(detections, p, wb, wbLo, wbHi)
%   [tracks, confirmedKeys] = trackLinkFwdBwd(detections, p, wb, wbLo, wbHi, pins)
%
%   pins : optional struct array of forced assignments (see trackLinkLAP),
%     anchored in ORIGINAL (forward-time) frame numbers. Applied directly
%     to the forward pass. For the backward pass (run on time-reversed
%     detections purely to build the confirmation map), each pin is
%     transformed into backward-frame coordinates so a manually-forced
%     link is confirmed by both passes rather than spuriously showing as
%     "unconfirmed" just because it was never discovered independently.
%
%   Runs trackLinkLAP twice — forward and backward in time — then classifies
%   each frame-to-frame link in the forward tracks as *confirmed* or
%   *unconfirmed* based on whether the backward pass made the same link.
%
%   Rationale
%   ---------
%   A wrong link (e.g. a fast mito grabbing a nearby slow one) tends to be
%   asymmetric: the forward pass makes it but the backward pass, starting
%   from the other end, takes a different path.  A correct link between the
%   same two objects will appear in both passes.  Confirmed links therefore
%   have higher confidence without any change to the cost model.
%
%   The backward pass is run on the time-reversed detection sequence.  A
%   backward link between label A at backward-frame k and label B at
%   backward-frame k+1 corresponds to a forward-time link between label B
%   at original frame nT+1-(k+1) and label A at original frame nT+1-k.
%   These pairs are stored in a containers.Map (confirmedKeys) for O(1)
%   lookup.
%
%   The `confirmed` field is NOT attached to tracks here — it is applied
%   after trackGapClose (which concatenates known fields by name and would
%   fail on unknown fields).  Call applyConfirmedField(tracks, confirmedKeys)
%   in trackCompute after gap closing.
%
%   INPUTS
%     detections      : [nT x 1] cell array from trackDetectionsExtract.
%     p               : parameter struct (see trackParamsDefault).
%     wb, wbLo, wbHi : optional waitbar handle and fractional range.
%                       Waitbar is split evenly between the two passes.
%
%   OUTPUTS
%     tracks        : standard struct array from trackLinkLAP (no extra fields).
%     confirmedKeys : containers.Map of confirmed link keys, used by
%                     applyConfirmedField.  Empty Map if no backward tracks.
%
%   See also: trackLinkLAP, trackCompute, applyConfirmedField

    nT    = numel(detections);
    useWb = nargin >= 3 && ~isempty(wb) && ishandle(wb);
    if nargin < 6, pins = []; end
    lo    = 0; hi = 1;
    if useWb, lo = wbLo; hi = wbHi; end
    mid   = lo + 0.5*(hi - lo);

    % ---- Forward pass -------------------------------------------------------
    if useWb
        waitbar(lo, wb, 'Forward linking...');
        tracks = trackLinkLAP(detections, p, wb, lo, mid, pins);
    else
        tracks = trackLinkLAP(detections, p, [], [], [], pins);
    end

    % ---- Backward pass on time-reversed detections --------------------------
    detBwd = flip(detections);   % detBwd{t} = detections{nT+1-t}

    % Transform forward-time pins into backward-frame coordinates. A pin
    % (frame=N, labelID=source@N, labelID2=target@N+1) corresponds, in
    % backward time, to a link starting at backward frame (nT-N) with the
    % TARGET label (the later original object becomes the backward
    % "source") and ending at backward frame (nT-N+1) with the SOURCE
    % label — i.e. source/target swap and the frame anchor shifts.
    pinsBwd = pins;
    for k = 1:numel(pins)
        pinsBwd(k).frame    = nT - pins(k).frame;
        pinsBwd(k).labelID  = pins(k).labelID2;
        pinsBwd(k).labelID2 = pins(k).labelID;
    end

    if useWb
        waitbar(mid, wb, 'Backward linking...');
        tracksBwd = trackLinkLAP(detBwd, p, wb, mid, hi, pinsBwd);
    else
        tracksBwd = trackLinkLAP(detBwd, p, [], [], [], pinsBwd);
    end

    % ---- Build confirmed-link key set from backward tracks ------------------
    %
    % A backward link at consecutive backward frames (fb_i, fb_i+1) between
    % labels (La, Lb) represents an original-time forward link:
    %
    %   Lb  at original frame  nT+1-fb_{i+1}  (earlier)
    %   La  at original frame  nT+1-fb_i       (later)
    %
    % Key format: '<labelEarly>_<frameEarly>_<labelLate>_<frameLate>'
    confirmedKeys = containers.Map('KeyType','char','ValueType','logical');

    for k = 1:numel(tracksBwd)
        tr = tracksBwd(k);
        for i = 1:numel(tr.frames)-1
            fbEarly = tr.frames(i+1);   % higher backward frame → earlier original
            fbLate  = tr.frames(i);     % lower  backward frame → later  original
            origEarly = nT + 1 - fbEarly;
            origLate  = nT + 1 - fbLate;
            labelEarly = tr.labelID(i+1);
            labelLate  = tr.labelID(i);
            key = linkKey(labelEarly, origEarly, labelLate, origLate);
            confirmedKeys(key) = true;
        end
    end
end

% -------------------------------------------------------------------------
function k = linkKey(labelA, frameA, labelB, frameB)
%LINKKEY  Canonical string key for a forward link A(frameA) → B(frameB).
    k = sprintf('%d_%d_%d_%d', labelA, frameA, labelB, frameB);
end
