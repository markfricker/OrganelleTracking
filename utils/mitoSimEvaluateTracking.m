function metrics = mitoSimEvaluateTracking(tracks, groundTruth, minTrackLength)
%MITOSIMALUATETRACKING  Compare tracking output against simulator ground truth.
%
%   metrics = mitoSimEvaluateTracking(tracks, groundTruth)
%   metrics = mitoSimEvaluateTracking(tracks, groundTruth, minTrackLength)
%
%   Uses nearest-neighbour matching between detected track centroids and the
%   ground-truth trajectories to compute per-track position error, recall,
%   and precision.
%
%   INPUTS
%     tracks         : struct array from trackCompute.
%     groundTruth    : cell array from mitoSimGenerate (one cell per frame,
%                      each containing a struct array with fields x, y, id).
%     minTrackLength : minimum number of frames a detected track must span
%                      to be counted (default: 5).
%
%   OUTPUT  metrics struct with fields:
%     nGT              number of ground-truth objects
%     nDetected        number of detected tracks (after length filter)
%     recall           fraction of GT objects matched (0..1)
%     precision        fraction of detected tracks matched to a GT object (0..1)
%     f1               2 * precision * recall / (precision + recall)
%     matchedIDs       [nMatched x 2] -- columns: [GT id, detected track index]
%     meanPosError     mean centroid error over matched frames (px)
%     medianPosError   median centroid error (px)
%     posErrorPerTrack [nMatched x 1] mean per-track position error (px)

    if nargin < 3 || isempty(minTrackLength)
        minTrackLength = 5;
    end

    nT  = numel(groundTruth);
    nGT = numel(groundTruth{1});  % objects present in first frame

    % Filter detected tracks by length.
    longTracks = tracks(arrayfun(@(t) numel(t.frames) >= minTrackLength, tracks));
    nDetected  = numel(longTracks);

    fprintf('\n=== Tracking Evaluation ===\n');
    fprintf('GT mitochondria:    %d\n', nGT);
    fprintf('Detected tracks:    %d  (>= %d frames)\n', nDetected, minTrackLength);

    if nDetected == 0 || nGT == 0
        metrics.nGT            = nGT;
        metrics.nDetected      = nDetected;
        metrics.recall         = 0;
        metrics.precision      = 0;
        metrics.f1             = 0;
        metrics.matchedIDs     = zeros(0,2);
        metrics.meanPosError   = NaN;
        metrics.medianPosError = NaN;
        metrics.posErrorPerTrack = [];
        return
    end

    % Build GT trajectories: cell{gtID} = [nT x 2] [x y] (NaN if absent).
    gtTraj = buildGTTrajectories(groundTruth, nGT, nT);

    % Match each detected track to the nearest GT trajectory by mean
    % nearest-frame displacement.
    matchThresh = 10;   % px -- maximum mean error to count as a true match
    assigned    = false(nGT, 1);
    matchedIDs  = zeros(0, 2);
    posErrors   = zeros(0, 1);

    for d = 1:nDetected
        tr    = longTracks(d);
        bestG = 0;
        bestE = Inf;

        for g = 1:nGT
            xy_gt  = gtTraj{g};
            frames = tr.frames;
            % Only compare frames where both exist.
            validF = frames(frames <= nT);
            if isempty(validF), continue, end
            gt_xy_f = xy_gt(validF, :);
            keep    = ~any(isnan(gt_xy_f), 2);
            if sum(keep) < 2, continue, end

            dx  = tr.x(keep)' - gt_xy_f(keep, 1);
            dy  = tr.y(keep)' - gt_xy_f(keep, 2);
            err = mean(hypot(dx, dy));

            if err < bestE
                bestE = err;
                bestG = g;
            end
        end

        if bestG > 0 && bestE < matchThresh && ~assigned(bestG)
            assigned(bestG) = true;
            matchedIDs(end+1, :) = [bestG, d]; %#ok<AGROW>
            posErrors(end+1)     = bestE;       %#ok<AGROW>
        end
    end

    nMatched          = size(matchedIDs, 1);
    metrics.nGT            = nGT;
    metrics.nDetected      = nDetected;
    metrics.recall         = nMatched / nGT;
    metrics.precision      = nMatched / max(1, nDetected);
    metrics.f1             = 2 * metrics.precision * metrics.recall / ...
                             max(eps, metrics.precision + metrics.recall);
    metrics.matchedIDs     = matchedIDs;
    metrics.meanPosError   = mean(posErrors);
    metrics.medianPosError = median(posErrors);
    metrics.posErrorPerTrack = posErrors(:);

    fprintf('Matched:            %d\n', nMatched);
    fprintf('Recall:             %.1f%%\n', metrics.recall   * 100);
    fprintf('Precision:          %.1f%%\n', metrics.precision * 100);
    fprintf('F1 score:           %.3f\n',   metrics.f1);
    if ~isempty(posErrors)
        fprintf('Mean pos error:     %.2f px\n', metrics.meanPosError);
        fprintf('Median pos error:   %.2f px\n', metrics.medianPosError);
    end
end

% -------------------------------------------------------------------------

function gtTraj = buildGTTrajectories(groundTruth, nGT, nT)
    % gtTraj{g} = [nT x 2] matrix; rows are [x y] for that GT object at each
    % frame; NaN where the object is not present.
    gtTraj = cell(nGT, 1);
    for g = 1:nGT
        gtTraj{g} = NaN(nT, 2);
    end
    for t = 1:nT
        mito = groundTruth{t};
        for i = 1:numel(mito)
            g = mito(i).id;
            if g <= nGT
                gtTraj{g}(t, :) = [mito(i).x, mito(i).y];
            end
        end
    end
end
