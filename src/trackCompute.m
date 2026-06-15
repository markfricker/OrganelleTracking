function [tracks, mergeInfo, labelStack] = trackCompute(labelStack, p)
%TRACKCOMPUTE  Top-level organelle tracking from a segmented label stack.
%   Third output labelStack is the (possibly merge-resolved) label stack —
%   identical to the input unless p.resolveMerges split fused blobs.
%
%   tracks = trackCompute(labelStack, p)
%
%   Pipeline (with progress waitbar):
%     1. trackDetectionsExtract  - regionprops features per frame  (0-40%)
%     2. trackLinkLAP            - frame-to-frame LAP linking       (40-90%)
%     3. trackGapClose           - second-pass global gap closing   (90-100%)
%
%   A waitbar with a Cancel button is shown for the duration.  Cancelling
%   at any stage returns an empty tracks array without throwing an error.
%
%   INPUTS
%     labelStack : [nY nX nT] uint16 label matrix.
%                 Each frame contains objects labelled 1..N (0 = background).
%                 Single channel, single Z assumed -- extract the desired
%                 C/Z plane before calling.
%     p          : parameter struct (see trackParamsDefault for defaults).
%                  Missing fields are filled from trackParamsDefault.
%
%   OUTPUT
%     tracks : struct array, one element per track (sorted by start frame).
%              Fields per track (all row vectors):
%
%       id            Unique integer.
%       frames        Frame indices (1-based).
%       x, y          Centroid (px).
%       vx, vy        Smoothed velocity (px/frame); NaN first observation.
%       speed         ||v|| in px/frame; NaN first observation.
%       speedUm       speed x pixelSize / frameInterval  (um/s).
%       area          Object area (px^2).
%       majorAxis     Major axis (px).
%       minorAxis     Minor axis (px).
%       aspectRatio   majorAxis / minorAxis.
%       orientation   Major-axis angle (deg, -90..90).
%       eccentricity  0 = circle, 1 = line.
%       labelID       Label value in labelStack (0 = gap-interpolated frame).
%       streaming     Logical flag: true at frames where speed > p.streamingThreshold.
%
%   Notes
%     Requires MATLAB Computer Vision Toolbox (assignDetectionsToTracks).

    p = mergeDefaults(p);

    nT = size(labelStack, 3);

    %% Create waitbar.
    wb = waitbar(0, sprintf('Tracking: extracting features (0 / %d frames)...', nT), ...
        'Name', 'Organelle Tracking', ...
        'CreateCancelBtn', 'setappdata(gcbf,''cancel'',true)');
    setappdata(wb, 'cancel', false);
    cleanupWb = onCleanup(@() closeWaitbar(wb));
    drawnow;       % force render before the computation loop
    figure(wb);    % bring waitbar to front (App Designer uifigure otherwise covers it)

    tracks = struct([]);
    mergeInfo = struct('nSplit', 0, 'frames', [], ...
        'splits', struct('frame', {}, 'hostLabel', {}, 'newLabel', {}, ...
                         'pHidden', {}, 'pHost', {}));

    try
        %% 1. Extract features (0 -> 0.40).
        detections = trackDetectionsExtract(labelStack, p, wb, 0.00, 0.40);

        %% 2. Frame-to-frame linking (0.40 -> 0.88).
        confirmedKeys = [];
        if getf(p, 'useFwdBwd', false)
            waitbar(0.40, wb, 'Linking detections (forward + backward)...');
            [tracks, confirmedKeys] = trackLinkFwdBwd(detections, p, wb, 0.40, 0.88);
        else
            waitbar(0.40, wb, 'Linking detections...');
            tracks = trackLinkLAP(detections, p, wb, 0.40, 0.88);
        end

        if isempty(tracks)
            waitbar(1, wb, 'Done -- no tracks found.');
            return
        end

        %% 3. Gap closing (0.88 -> 0.96).
        waitbar(0.88, wb, sprintf('Closing gaps across %d segments...', numel(tracks)));
        tracks = trackGapClose(tracks, detections, p);

        if isempty(tracks)
            waitbar(1, wb, 'Done -- no tracks survived filtering.');
            return
        end

        %% 4. Apply confirmed field from fwd-bwd pass (if used).
        if ~isempty(confirmedKeys)
            waitbar(0.96, wb, 'Applying confirmation flags...');
            tracks = applyConfirmedField(tracks, confirmedKeys);
        end

        %% 4b. Resolve transient merges (optional): split fused blobs using
        %% temporally-interpolated seeds, then re-track on the updated stack.
        if p.resolveMerges
            waitbar(0.96, wb, 'Resolving transient merges...');
            [labelStack, mInfo] = trackResolveMerges(labelStack, tracks, p);
            mergeInfo = mInfo;
            if mInfo.nSplit > 0
                detections = trackDetectionsExtract(labelStack, p, wb, 0.96, 0.97);
                if p.useFwdBwd
                    [tracks, confirmedKeys] = trackLinkFwdBwd(detections, p, wb, 0.97, 0.985);
                else
                    tracks = trackLinkLAP(detections, p, wb, 0.97, 0.985);
                    confirmedKeys = [];
                end
                tracks = trackGapClose(tracks, detections, p);
                if ~isempty(confirmedKeys)
                    tracks = applyConfirmedField(tracks, confirmedKeys);
                end
            end
        end

        %% 5. Derived quantities (0.96 -> 1.00).
        waitbar(0.96, wb, 'Computing velocities and streaming flags...');
        umPerPxPerFrame = p.pixelSize / p.frameInterval;
        for i = 1:numel(tracks)
            tracks(i).speedUm   = tracks(i).speed * umPerPxPerFrame;
            tracks(i).streaming = tracks(i).speed > p.streamingThreshold;
        end

        %% 6. Sort by start frame.
        [~, ord] = sort(arrayfun(@(t) t.frames(1), tracks));
        tracks = tracks(ord);

        waitbar(1, wb, sprintf('Done -- %d tracks found.', numel(tracks)));
        pause(0.3);   % brief pause so user sees the final message

    catch ME
        if strcmp(ME.identifier, 'trackCompute:cancelled')
            % User pressed Cancel -- return empty silently.
            tracks = struct([]);
        else
            rethrow(ME);
        end
    end
end

% -------------------------------------------------------------------------

function closeWaitbar(wb)
    try
        if ishandle(wb)
            delete(wb);
        end
    catch
    end
end

function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field)), v = s.(field); else, v = default; end
end

function p = mergeDefaults(p)
    d = trackParamsDefault();
    fns = fieldnames(d);
    for k = 1:numel(fns)
        if ~isfield(p, fns{k})
            p.(fns{k}) = d.(fns{k});
        end
    end
    if ~isfield(p, 'streamingThreshold')
        p.streamingThreshold = 11;   % px/frame (~1 um/s at 90 nm/px, 0.4 s/frame)
    end
end
