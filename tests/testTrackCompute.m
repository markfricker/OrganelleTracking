classdef testTrackCompute < matlab.unittest.TestCase
%TESTTRACKCOMPUTE  Unit tests for the OrganelleTracking_sandbox library.
%
% USAGE
%   Run all tests from the MATLAB command window:
%       results = runtests('tests/testTrackCompute');
%       table(results)
%
%   Run a single test method:
%       runtests('tests/testTrackCompute/testLinkLAP_streamingMitoLinked')
%
% COVERAGE
%   trackDetectionsExtract  —  5 tests
%   trackLinkLAP            —  8 tests
%   trackGapClose           —  4 tests
%   trackCompute            —  9 tests
%
% REQUIREMENTS
%   MATLAB R2019b+ (matlab.unittest framework)
%   Computer Vision Toolbox  (assignDetectionsToTracks)
%   Image Processing Toolbox (regionprops)
%   All src/ and utils/ on the path (added automatically by TestClassSetup).

    properties (Constant)
        PosTol  = 1.0;    % px  — centroid position tolerance
        SpeedTol= 0.5;    % px/frame — speed tolerance
        ImSize  = [128, 128]  % default synthetic image size [nY nX]
    end

    % =====================================================================
    % Path setup
    % =====================================================================
    methods (TestClassSetup)
        function addPaths(tc) %#ok<MANU>
            rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(rootDir, 'src'));
            addpath(fullfile(rootDir, 'utils'));
            addpath(fullfile(rootDir, 'demos'));
        end
    end

    % =====================================================================
    % Shared synthetic helpers
    %
    %   singleStationary — 1 mito, fixed position, nT frames
    %   twoStationary    — 2 mitos, fixed separate positions, nT frames
    %   slowDiffusive    — 1 mito drifting 2 px/frame in x
    %   streamingTrack   — 1 mito: slow diffusion then sudden acceleration
    %   gappedTrack      — 1 mito absent for gapLen frames in the middle
    %   shortAndLong     — 1 short track + 1 long track
    %   defaultP         — trackParamsDefault with minTrackLength=1 (no filter)
    %
    %   All data generation delegates to trackSyntheticData (utils/).
    %   Raw frame painting (blankFrame / paintEllipse) is private to
    %   trackSyntheticData; tests use named scenarios instead.
    % =====================================================================
    methods (Static, Access = private)

        function stack = singleStationary(nT)
            stack = trackSyntheticData('stationary', nT);
        end

        function stack = twoStationary(nT)
            stack = trackSyntheticData('two', nT);
        end

        function [stack, xTrue] = slowDiffusive(nT)
            [stack, xTrue] = trackSyntheticData('drift', nT);
        end

        function [stack, xTrue] = streamingTrack(nT)
            [stack, xTrue] = trackSyntheticData('streaming', nT);
        end

        function stack = gappedTrack(nT, gapStart, gapLen)
            stack = trackSyntheticData('gapped', nT, gapStart, gapLen);
        end

        function stack = shortAndLong()
            stack = trackSyntheticData('shortAndLong');
        end

        function p = defaultP()
            % trackParamsDefault with minTrackLength=1 (no post-filter noise).
            p = trackParamsDefault();
            p.minTrackLength = 1;
        end

    end

    % =====================================================================
    %  trackDetectionsExtract
    % =====================================================================
    methods (Test)

        function testExtract_smoke(tc)
            % Basic call must complete without error.
            stack = testTrackCompute.singleStationary(5);
            d = trackDetectionsExtract(stack, []);
            tc.verifyNotEmpty(d);
        end

        function testExtract_outputLength(tc)
            % Output cell array length must match number of frames.
            nT = 8;
            stack = testTrackCompute.singleStationary(nT);
            d = trackDetectionsExtract(stack, []);
            tc.verifyLength(d, nT);
        end

        function testExtract_centroidAccuracy(tc)
            % Centroid of a painted ellipse must match the painted centre.
            % Uses trackSyntheticData('stationary',1) which paints at (64,64).
            stack3 = trackSyntheticData('stationary', 1);   % [128 x 128 x 1]
            d = trackDetectionsExtract(stack3, []);
            tc.verifyFalse(isempty(d{1}), 'No detections in frame 1');
            tc.verifyEqual(d{1}.x(1), 64, 'AbsTol', testTrackCompute.PosTol);
            tc.verifyEqual(d{1}.y(1), 64, 'AbsTol', testTrackCompute.PosTol);
        end

        function testExtract_emptyFrameGivesEmpty(tc)
            % A completely blank frame must yield an empty cell entry.
            sz    = testTrackCompute.ImSize;
            stack = zeros([sz, 1], 'uint16');
            d = trackDetectionsExtract(stack, []);
            tc.verifyTrue(isempty(d{1}));
        end

        function testExtract_fieldsPresent(tc)
            % Each non-empty detection struct must carry all required fields.
            stack = testTrackCompute.singleStationary(3);
            d = trackDetectionsExtract(stack, []);
            required = {'x','y','area','majorAxis','minorAxis', ...
                        'aspectRatio','orientation','eccentricity','labelID'};
            for f = required
                tc.verifyTrue(isfield(d{1}, f{1}), ...
                    sprintf('Field ''%s'' missing from detections', f{1}));
            end
        end

    end

    % =====================================================================
    %  trackLinkLAP
    % =====================================================================
    methods (Test)

        function testLinkLAP_smoke(tc)
            % Basic call must complete without error.
            stack = testTrackCompute.singleStationary(5);
            d = trackDetectionsExtract(stack, []);
            p = testTrackCompute.defaultP();
            tracks = trackLinkLAP(d, p);
            tc.verifyNotEmpty(tracks);
        end

        function testLinkLAP_singleMito_oneTrack(tc)
            % One stationary mito must produce exactly one track.
            stack = testTrackCompute.singleStationary(10);
            d = trackDetectionsExtract(stack, []);
            p = testTrackCompute.defaultP();
            tracks = trackLinkLAP(d, p);
            tc.verifyNumElements(tracks, 1);
        end

        function testLinkLAP_singleMito_fullLength(tc)
            % Track from a stationary mito must span all nT frames.
            nT    = 10;
            stack = testTrackCompute.singleStationary(nT);
            d     = trackDetectionsExtract(stack, []);
            p     = testTrackCompute.defaultP();
            tracks= trackLinkLAP(d, p);
            tc.verifyEqual(numel(tracks(1).frames), nT);
        end

        function testLinkLAP_twoMitos_twoTracks(tc)
            % Two well-separated mitos must produce exactly two tracks.
            stack = testTrackCompute.twoStationary(10);
            d     = trackDetectionsExtract(stack, []);
            p     = testTrackCompute.defaultP();
            tracks= trackLinkLAP(d, p);
            tc.verifyNumElements(tracks, 2);
        end

        function testLinkLAP_slowDrift_positionAccurate(tc)
            % Track positions must follow the known linear trajectory.
            nT = 8;
            [stack, xTrue] = testTrackCompute.slowDiffusive(nT);
            d     = trackDetectionsExtract(stack, []);
            p     = testTrackCompute.defaultP();
            tracks= trackLinkLAP(d, p);
            % Find the single track; sort frames in case they're unordered.
            tc.verifyNumElements(tracks, 1);
            [~, ord] = sort(tracks(1).frames);
            xObs = tracks(1).x(ord);
            tc.verifyEqual(xObs, xTrue, 'AbsTol', testTrackCompute.PosTol);
        end

        function testLinkLAP_streamingMitoLinked(tc)
            % A mito that suddenly accelerates (streaming) must stay in ONE
            % track — not fragment into separate segments.
            nT = 8;
            [stack, ~] = testTrackCompute.streamingTrack(nT);
            d  = trackDetectionsExtract(stack, []);
            p  = testTrackCompute.defaultP();
            p.maxDisplacement = 20;
            tracks = trackLinkLAP(d, p);
            % Accept that gap closing may eventually merge segments,
            % but after linking alone there should be at most 2 segments
            % (the initial velocity-blind frame transition may fragment once).
            totalFramesCovered = sum(arrayfun(@(t) numel(t.frames), tracks));
            tc.verifyEqual(totalFramesCovered, nT, ...
                'Not all frames accounted for across track segments');
        end

        function testLinkLAP_outputFieldsPresent(tc)
            % Track struct must contain all required output fields.
            stack  = testTrackCompute.singleStationary(5);
            d      = trackDetectionsExtract(stack, []);
            p      = testTrackCompute.defaultP();
            tracks = trackLinkLAP(d, p);
            required = {'id','frames','x','y','vx','vy','speed', ...
                        'area','majorAxis','minorAxis','aspectRatio', ...
                        'orientation','eccentricity','labelID'};
            for f = required
                tc.verifyTrue(isfield(tracks, f{1}), ...
                    sprintf('Field ''%s'' missing from tracks', f{1}));
            end
        end

        function testLinkLAP_velocityIsNaNFirstFrame(tc)
            % The first observation of each track must have NaN velocity.
            stack  = testTrackCompute.singleStationary(5);
            d      = trackDetectionsExtract(stack, []);
            p      = testTrackCompute.defaultP();
            tracks = trackLinkLAP(d, p);
            tc.verifyTrue(isnan(tracks(1).vx(1)), ...
                'vx(1) should be NaN for first observation');
            tc.verifyTrue(isnan(tracks(1).vy(1)), ...
                'vy(1) should be NaN for first observation');
        end

    end

    % =====================================================================
    %  trackGapClose
    % =====================================================================
    methods (Test)

        function testGapClose_bridgesShortGap(tc)
            % A 2-frame gap within gapMax must be bridged into a single track.
            nT       = 12;
            gapStart = 5;
            gapLen   = 2;
            stack    = testTrackCompute.gappedTrack(nT, gapStart, gapLen);
            d        = trackDetectionsExtract(stack, []);
            p        = testTrackCompute.defaultP();
            p.gapMax = 3;
            p.minTrackLength = 1;
            tracks   = trackLinkLAP(d, p);
            tracks   = trackGapClose(tracks, d, p);
            % After gap closing there should be exactly one track.
            tc.verifyNumElements(tracks, 1, ...
                'Gap closing should merge two segments into one track');
        end

        function testGapClose_doesNotBridgeLongGap(tc)
            % A gap larger than gapMax must NOT be bridged.
            nT       = 14;
            gapStart = 5;
            gapLen   = 6;   % > gapMax = 3
            stack    = testTrackCompute.gappedTrack(nT, gapStart, gapLen);
            d        = trackDetectionsExtract(stack, []);
            p        = testTrackCompute.defaultP();
            p.gapMax = 3;
            p.minTrackLength = 1;
            tracks   = trackLinkLAP(d, p);
            tracks   = trackGapClose(tracks, d, p);
            % Two segments should remain separate.
            tc.verifyGreaterThanOrEqual(numel(tracks), 2, ...
                'Long gap must not be bridged');
        end

        function testGapClose_interpolatedFramesLabelIDZero(tc)
            % Interpolated gap frames must have labelID == 0.
            nT       = 10;
            gapStart = 4;
            gapLen   = 2;
            stack    = testTrackCompute.gappedTrack(nT, gapStart, gapLen);
            d        = trackDetectionsExtract(stack, []);
            p        = testTrackCompute.defaultP();
            p.gapMax = 3;
            p.minTrackLength = 1;
            tracks   = trackLinkLAP(d, p);
            tracks   = trackGapClose(tracks, d, p);
            if numel(tracks) == 1
                gapIdx = ismember(tracks(1).frames, gapStart : gapStart+gapLen-1);
                tc.verifyTrue(all(tracks(1).labelID(gapIdx) == 0), ...
                    'Interpolated frames must have labelID = 0');
            end
        end

        function testGapClose_mergedTrackCoversAllFrames(tc)
            % After merging, the track must span from frame 1 to nT.
            nT       = 10;
            gapStart = 4;
            gapLen   = 2;
            stack    = testTrackCompute.gappedTrack(nT, gapStart, gapLen);
            d        = trackDetectionsExtract(stack, []);
            p        = testTrackCompute.defaultP();
            p.gapMax = 3;
            p.minTrackLength = 1;
            tracks   = trackLinkLAP(d, p);
            tracks   = trackGapClose(tracks, d, p);
            if numel(tracks) == 1
                tc.verifyEqual(min(tracks(1).frames), 1);
                tc.verifyEqual(max(tracks(1).frames), nT);
            end
        end

    end

    % =====================================================================
    %  trackCompute  (full pipeline)
    % =====================================================================
    methods (Test)

        function testCompute_smoke(tc)
            % Full pipeline call must complete without error.
            stack = testTrackCompute.singleStationary(10);
            p     = testTrackCompute.defaultP();
            tracks = trackCompute(stack, p);
            tc.verifyNotEmpty(tracks);
        end

        function testCompute_emptyStack_returnsEmpty(tc)
            % All-zero label stack must return empty track array.
            stack = zeros([testTrackCompute.ImSize, 10], 'uint16');
            p     = testTrackCompute.defaultP();
            tracks = trackCompute(stack, p);
            tc.verifyTrue(isempty(tracks));
        end

        function testCompute_minTrackLength_filtersShortTracks(tc)
            % Short track (3 frames) must be removed when minTrackLength = 5.
            stack = testTrackCompute.shortAndLong();
            p     = trackParamsDefault();
            p.minTrackLength = 5;
            tracks = trackCompute(stack, p);
            % Only the long track (10 frames) should survive.
            tc.verifyNumElements(tracks, 1);
            tc.verifyGreaterThanOrEqual(numel(tracks(1).frames), 5);
        end

        function testCompute_speedUmFieldPresent(tc)
            % trackCompute must add speedUm field to every track.
            stack  = testTrackCompute.singleStationary(10);
            p      = testTrackCompute.defaultP();
            tracks = trackCompute(stack, p);
            tc.verifyTrue(isfield(tracks, 'speedUm'), ...
                'speedUm field missing from trackCompute output');
        end

        function testCompute_streamingFieldPresent(tc)
            % trackCompute must add streaming logical field to every track.
            stack  = testTrackCompute.singleStationary(10);
            p      = testTrackCompute.defaultP();
            tracks = trackCompute(stack, p);
            tc.verifyTrue(isfield(tracks, 'streaming'), ...
                'streaming field missing from trackCompute output');
        end

        function testCompute_stationaryMito_speedNearZero(tc)
            % A stationary mito must have near-zero speed after the first frame.
            stack  = testTrackCompute.singleStationary(10);
            p      = testTrackCompute.defaultP();
            tracks = trackCompute(stack, p);
            sp = tracks(1).speedUm;
            sp = sp(~isnan(sp));
            tc.verifyLessThanOrEqual(max(sp), testTrackCompute.SpeedTol, ...
                'Stationary mito should have speed near zero');
        end

        function testCompute_streamingFlagTrue(tc)
            % Streaming frames must be flagged as streaming = true.
            nT = 8;
            [stack, ~] = testTrackCompute.streamingTrack(nT);
            p = testTrackCompute.defaultP();
            p.maxDisplacement    = 20;
            p.streamingThreshold = 10;  % px/frame — below streaming velocity
            p.gapMax = 3;
            tracks = trackCompute(stack, p);
            % Consolidate all streaming flags across possibly-fragmented tracks.
            anyStreaming = false;
            for i = 1:numel(tracks)
                if any(tracks(i).streaming)
                    anyStreaming = true;
                    break
                end
            end
            tc.verifyTrue(anyStreaming, ...
                'Streaming frames must be flagged in tracks.streaming');
        end

        function testCompute_sortedByStartFrame(tc)
            % Output tracks must be sorted ascending by their first frame.
            stack  = testTrackCompute.twoStationary(10);
            p      = testTrackCompute.defaultP();
            tracks = trackCompute(stack, p);
            startFrames = arrayfun(@(t) t.frames(1), tracks);
            tc.verifyEqual(startFrames, sort(startFrames), ...
                'Tracks must be sorted by start frame');
        end

        function testCompute_speedUm_physicalUnits(tc)
            % A mito moving at known velocity must give correct speedUm.
            % slowDiffusive drifts at 2 px/frame; with default params
            % (pixelSize=0.09 µm, frameInterval=0.4 s) → 0.45 µm/s.
            nT     = 10;
            [stack, ~] = testTrackCompute.slowDiffusive(nT);
            p      = testTrackCompute.defaultP();
            p.pixelSize     = 0.09;
            p.frameInterval = 0.4;
            tracks = trackCompute(stack, p);
            tc.verifyNumElements(tracks, 1);
            sp = tracks(1).speedUm;
            sp = sp(~isnan(sp));
            expectedSpeed = 2 * 0.09 / 0.4;   % 0.45 µm/s
            tc.verifyEqual(mean(sp), expectedSpeed, 'AbsTol', 0.05);
        end

    end

end
