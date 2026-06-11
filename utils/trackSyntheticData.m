function varargout = trackSyntheticData(scenario, varargin)
%TRACKSYNTHETICDATA  Generate synthetic label stacks for tracking demos/tests.
%
%   stack                  = trackSyntheticData('stationary',  nT)
%   stack                  = trackSyntheticData('two',         nT)
%   [stack, xTrue]         = trackSyntheticData('drift',       nT)
%   [stack, xTrue]         = trackSyntheticData('streaming',   nT)
%   stack                  = trackSyntheticData('gapped',      nT, gapStart, gapLen)
%   stack                  = trackSyntheticData('shortAndLong')
%
%   All stacks are [128 x 128 x nT] uint16.  Objects are 8x4 px filled
%   ellipses (major axis horizontal, orientation = 0 deg).
%
%   Scenarios
%   ---------
%   'stationary'    One mito fixed at (64,64) for nT frames.
%   'two'           Two mitos at (40,64) and (80,64), both present all frames.
%   'drift'         One mito drifting 2 px/frame to the right from x=30.
%                   Second output xTrue gives the ground-truth x positions.
%   'streaming'     One mito: slow drift (2 px/frame) for frames 1-4, then
%                   sudden streaming (15 px/frame) from frame 5 onward.
%                   Second output xTrue gives ground-truth x positions.
%   'gapped'        One mito at (64,64) absent for gapLen frames starting at
%                   gapStart.  Requires extra arguments: nT, gapStart, gapLen.
%   'shortAndLong'  Short track (frames 1-3) at (30,30) plus long track
%                   (frames 1-10) at (90,90).  No nT argument needed.
%
%   EXAMPLE
%     stack = trackSyntheticData('streaming', 30);
%     p = trackParamsDefault();
%     p.minTrackLength = 1;
%     tracks = trackCompute(stack, p);
%     trackMovie(stack, [], tracks, p);

    imSize = [128, 128];

    switch lower(scenario)

        case 'stationary'
            nT    = getArg(varargin, 1, 20);
            stack = makeStack(imSize, nT);
            for t = 1:nT
                stack(:,:,t) = paintEllipse(stack(:,:,t), 64, 64, 8, 4, 0, 1);
            end
            varargout = {stack};

        case 'two'
            nT    = getArg(varargin, 1, 20);
            stack = makeStack(imSize, nT);
            for t = 1:nT
                f = blankFrame(imSize);
                f = paintEllipse(f, 40, 64, 8, 4, 0, 1);
                f = paintEllipse(f, 80, 64, 8, 4, 0, 2);
                stack(:,:,t) = f;
            end
            varargout = {stack};

        case 'drift'
            nT    = getArg(varargin, 1, 20);
            stack = makeStack(imSize, nT);
            xTrue = zeros(1, nT);
            cx = 30;
            for t = 1:nT
                cx = cx + 2;
                xTrue(t) = cx;
                stack(:,:,t) = paintEllipse(blankFrame(imSize), cx, 64, 8, 4, 0, 1);
            end
            varargout = {stack, xTrue};

        case 'streaming'
            nT    = getArg(varargin, 1, 20);
            stack = makeStack(imSize, nT);
            xTrue = zeros(1, nT);
            cx = 30;
            for t = 1:nT
                if t <= 4
                    cx = cx + 2;
                else
                    cx = cx + 15;
                end
                xTrue(t) = cx;
                stack(:,:,t) = paintEllipse(blankFrame(imSize), cx, 64, 8, 4, 0, 1);
            end
            varargout = {stack, xTrue};

        case 'gapped'
            nT       = getArg(varargin, 1, 20);
            gapStart = getArg(varargin, 2, 6);
            gapLen   = getArg(varargin, 3, 2);
            stack    = makeStack(imSize, nT);
            for t = 1:nT
                if t >= gapStart && t < gapStart + gapLen
                    continue
                end
                stack(:,:,t) = paintEllipse(blankFrame(imSize), 64, 64, 8, 4, 0, 1);
            end
            varargout = {stack};

        case 'shortandlong'
            nT    = 10;
            stack = makeStack(imSize, nT);
            for t = 1:nT
                f = blankFrame(imSize);
                if t <= 3
                    f = paintEllipse(f, 30, 30, 8, 4, 0, 1);
                end
                f = paintEllipse(f, 90, 90, 8, 4, 0, 2);
                stack(:,:,t) = f;
            end
            varargout = {stack};

        otherwise
            error('trackSyntheticData:unknownScenario', ...
                'Unknown scenario ''%s''. Choose: stationary, two, drift, streaming, gapped, shortAndLong.', ...
                scenario);
    end
end

% -------------------------------------------------------------------------

function stack = makeStack(imSize, nT)
    stack = zeros([imSize, nT], 'uint16');
end

function frame = blankFrame(imSize)
    frame = zeros(imSize(1), imSize(2), 'uint16');
end

function frame = paintEllipse(frame, cx, cy, a, b, orientDeg, labelVal)
    [nY, nX] = size(frame);
    [xx, yy] = meshgrid(1:nX, 1:nY);
    theta  = deg2rad(orientDeg);
    xr =  (xx - cx) .* cos(theta) + (yy - cy) .* sin(theta);
    yr = -(xx - cx) .* sin(theta) + (yy - cy) .* cos(theta);
    inside = (xr / a).^2 + (yr / b).^2 <= 1;
    frame(inside) = uint16(labelVal);
end

function val = getArg(args, idx, default)
    if numel(args) >= idx && ~isempty(args{idx})
        val = args{idx};
    else
        val = default;
    end
end
