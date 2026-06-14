function tracks = applyConfirmedField(tracks, confirmedKeys)
%APPLYCONFIRMEDFIELD  Add 'confirmed' field to tracks from a fwd-bwd link key set.
%
%   tracks = applyConfirmedField(tracks, confirmedKeys)
%
%   Called in trackCompute after trackGapClose, so gap-interpolated frames
%   (labelID == 0) are already present.  Those frames are always marked
%   confirmed = false regardless of the key set.
%
%   For each direct detection frame (labelID > 0), the link FROM the
%   previous observed frame TO this frame is looked up in confirmedKeys.
%   If found, confirmed = true.  The first observation of every track has
%   no preceding link and is marked confirmed = true (it was detected, not
%   inferred).
%
%   INPUTS
%     tracks        : struct array from trackGapClose.
%     confirmedKeys : containers.Map returned by trackLinkFwdBwd.
%
%   OUTPUT
%     tracks : same struct array with a 'confirmed' logical row-vector field
%              added to every element.

    for k = 1:numel(tracks)
        tr   = tracks(k);
        nObs = numel(tr.frames);
        conf = false(1, nObs);

        % Find indices of real detections (not gap-interpolated).
        realIdx = find(tr.labelID > 0);

        % First real detection: confirmed by definition (it was detected).
        if ~isempty(realIdx)
            conf(realIdx(1)) = true;
        end

        % Subsequent real detections: look up the link from the previous
        % real detection in confirmedKeys.
        for r = 2:numel(realIdx)
            iNow  = realIdx(r);
            iPrev = realIdx(r-1);
            key = sprintf('%d_%d_%d_%d', ...
                tr.labelID(iPrev), tr.frames(iPrev), ...
                tr.labelID(iNow),  tr.frames(iNow));
            if isKey(confirmedKeys, key)
                conf(iNow) = true;
            end
        end

        tracks(k).confirmed = conf;
    end
end
