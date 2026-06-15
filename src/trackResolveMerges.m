function [labelStack, info] = trackResolveMerges(labelStack, tracks, p)
%TRACKRESOLVEMERGES  Split transiently-fused objects using temporal seeds.
%
%   [labelStack, info] = trackResolveMerges(labelStack, tracks, p)
%
%   Resolves "transient merges": frames where two mitochondria are segmented
%   as a single blob but are separately detected in flanking frames.
%
%   KEY IDEA
%   --------
%   After gap closing, the absorbed object's track is bridged across the
%   fused interval and its centroid is LINEARLY INTERPOLATED through those
%   frames (those observations carry labelID == 0).  Such an interpolated
%   centroid that falls INSIDE another object's blob is the hidden object's
%   temporal seed.  The blob is then split by a mask-only distance watershed
%   seeded at (a) the hidden object's interpolated position and (b) the
%   geodesic-farthest interior point (a proxy for the host object), and the
%   hidden portion is relabelled.  No intensity image is required.
%
%   The updated labelStack should be RE-TRACKED by the caller so the new
%   labels are linked into proper tracks (trackCompute does this when
%   p.resolveMerges is true).
%
%   INPUTS
%     labelStack : [nY nX nT] uint16 label matrix (gap-closed tracking input).
%     tracks     : struct array from trackGapClose (must contain labelID==0
%                  rows for gap-interpolated frames).
%     p          : parameter struct.  Fields used:
%                    mergeMinSeedSep  min geodesic separation (px) between the
%                                     two seeds for a split to be attempted
%                                     (default 4) — guards against splitting
%                                     small/compact blobs.
%
%   OUTPUTS
%     labelStack : updated stack with transiently-merged blobs split; each
%                  hidden object receives a fresh label in its frame.
%     info       : struct with fields
%                    nSplit  number of blob splits performed.
%                    frames  frame indices where a split occurred.
%
%   LIMITATIONS (v1)
%     - Only TRANSIENT merges (gap-bridgeable: fusion duration <= p.gapMax,
%       plausible motion).  Persistent fusion is untouched (no interpolated
%       seed).  This is intended.
%     - Two-way merges; 3+-way clumps are only partially handled.
%     - The host seed is a geodesic-farthest-point proxy, not the host track's
%       extrapolated position; adequate for adjacent bodies.  A bwd-identity
%       host estimate is a planned v2 refinement.
%
%   See also: trackGapClose, trackCompute, watershed, bwdistgeodesic

    [nY, nX, ~] = size(labelStack);
    info.nSplit = 0;
    info.frames = [];
    info.splits = struct('frame', {}, 'hostLabel', {}, 'newLabel', {}, ...
                         'pHidden', {}, 'pHost', {});

    minSep = getf(p, 'mergeMinSeedSep', 4);

    for t = 1:numel(tracks)
        tr = tracks(t);
        for i = 1:numel(tr.frames)

            % Trigger: a gap-interpolated observation (no real detection).
            if tr.labelID(i) ~= 0
                continue
            end

            f  = tr.frames(i);
            cH = round(tr.x(i));      % hidden-object seed, x = column
            rH = round(tr.y(i));      %                     y = row
            if rH < 1 || cH < 1 || rH > nY || cH > nX
                continue
            end

            frame = labelStack(:, :, f);
            L = frame(rH, cH);
            if L == 0
                continue        % interpolated point in empty space → plain gap, not a merge
            end

            blob = (frame == L);

            % Host seed = interior point geodesically farthest from the hidden
            % seed (opposite end of the fused blob).
            Dgeo = bwdistgeodesic(blob, cH, rH, 'quasi-euclidean');
            Dgeo(~blob) = -Inf;
            [farVal, farLin] = max(Dgeo(:));
            if ~isfinite(farVal) || farVal < minSep
                continue        % blob too small/compact to be two objects
            end
            [rHost, cHost] = ind2sub([nY, nX], farLin);

            % Mask-only split: impose the two seeds as minima of -distance.
            Dd = bwdist(~blob);
            seeds = false(nY, nX);
            seeds(rH, cH)       = true;
            seeds(rHost, cHost) = true;
            relief = imimposemin(-Dd, seeds);
            parts  = watershed(relief);
            parts(~blob) = 0;

            labHidden = parts(rH, cH);
            if labHidden == 0
                continue        % seed landed on a watershed ridge; skip
            end

            % Relabel the hidden portion with a fresh label; host keeps L.
            newLabel = max(frame(:)) + 1;
            frame(blob & parts == labHidden) = newLabel;
            labelStack(:, :, f) = frame;

            info.nSplit = info.nSplit + 1;
            info.frames(end+1) = f;
            info.splits(end+1) = struct('frame', f, 'hostLabel', double(L), ...
                'newLabel', double(newLabel), 'pHidden', [cH rH], 'pHost', [cHost rHost]);
        end
    end
end

% -------------------------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
