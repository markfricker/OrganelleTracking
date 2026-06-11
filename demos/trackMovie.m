function trackMovie(labelStack, imRaw, tracks, p, varargin)
%TRACKMOVIE  Animated frame-by-frame overlay of tracking results.
%
%   trackMovie(labelStack, imRaw, tracks, p)
%   trackMovie(labelStack, imRaw, tracks, p, Name, Value, ...)
%
%   Each frame shows:
%     - Raw image (grey, contrast-stretched once from full stack)
%     - Mitochondrion outlines from the label stack (white, current frame)
%     - Track trails: last TrailLength frames drawn as coloured lines;
%       streaming segments shown in red regardless of track colour.
%     - Current position: filled circle at track centroid.
%     - Track ID labels (optional).
%
%   Name-Value options
%   ------------------
%   TrailLength   Frames of history shown as a fading trail.  (default 15)
%   ColorBy       'track'  - one fixed colour per track (default)
%                 'speed'  - colour encodes instantaneous speed (blue-red)
%   ShowIDs       true/false - show track ID text.             (default true)
%   FrameRate     Display/export frame rate (fps).             (default 10)
%   OutputFile    Path to save video (.avi or .mp4).
%                 If empty or omitted, animates in a figure only.
%   FramePause    Pause between frames when animating (s).     (default 0.05)
%                 Set to 0 for maximum speed (limited by render time).
%
%   INPUTS
%     labelStack : [nY nX nT] uint16 label matrix from segmentation.
%     imRaw      : [nY nX nT] numeric raw intensity (any class).
%                  Pass [] or labelStack to fall back to label display.
%     tracks     : struct array from trackCompute.
%     p          : parameter struct used for tracking (for streamingThreshold).
%                  Pass [] to use trackParamsDefault.
%
%   EXAMPLE
%     p = trackParamsDefault();
%     tracks = trackCompute(myLabels, p);
%     trackMovie(myLabels, myRaw, tracks, p, 'OutputFile', 'mito_tracks.mp4');

    %% Parse inputs.
    opts = parseOptions(varargin);

    if isempty(p)
        p = trackParamsDefault();
    end
    if ~isfield(p, 'streamingThreshold')
        p.streamingThreshold = 11;
    end

    [nY, nX, nT] = size(labelStack);

    % Fall back to label image if no raw provided.
    if isempty(imRaw)
        imRaw = single(labelStack > 0);
    end

    %% Pre-compute background: per-frame normalised intensity [nY nX nT].
    % Stretch contrast once from the full stack so brightness is consistent.
    imF    = double(imRaw);
    imLow  = prctile(imF(:), 1);
    imHigh = prctile(imF(:), 99);
    imF    = (imF - imLow) / max(imHigh - imLow, eps);
    imF    = max(0, min(1, imF));   % clamp

    %% Pre-compute object outlines per frame (for speed).
    outlines = precomputeOutlines(labelStack, nT);

    %% Colour assignment.
    nTracks = numel(tracks);
    trackColors = lines(max(nTracks, 1));
    trackColors = trackColors(mod((1:nTracks)-1, size(trackColors,1))+1, :);

    % Speed colormap (blue -> cyan -> yellow -> red) for ColorBy='speed'.
    if strcmp(opts.colorBy, 'speed')
        allSpeeds = [];
        for i = 1:nTracks
            sp = tracks(i).speedUm;
            allSpeeds = [allSpeeds, sp(~isnan(sp))]; %#ok<AGROW>
        end
        speedMax = prctile(allSpeeds, 98);
        if speedMax == 0, speedMax = 1; end
    else
        speedMax = 1;   % unused
    end

    %% Build track lookup: for each track, map frame -> (x, y, speed, streaming).
    % Store as cell array indexed by track number.
    trackFrameMap = buildFrameMap(tracks, nT);

    %% Set up figure.
    fig = figure('Name', 'Tracking movie', ...
                 'NumberTitle', 'off', ...
                 'Color', 'k', ...
                 'Units', 'pixels');
    % Size figure to match image aspect ratio (max 900 px wide).
    scale    = min(900 / nX, 700 / nY);
    figW     = round(nX * scale);
    figH     = round(nY * scale);
    fig.Position(3:4) = [figW, figH];

    ax = axes('Parent', fig, 'Position', [0 0 1 1]);
    axis(ax, 'off');

    %% Set up video writer if requested.
    writer = [];
    if ~isempty(opts.outputFile)
        [~, ~, ext] = fileparts(opts.outputFile);
        if strcmpi(ext, '.mp4')
            writer = VideoWriter(opts.outputFile, 'MPEG-4');
        else
            writer = VideoWriter(opts.outputFile, 'Motion JPEG AVI');
        end
        writer.FrameRate = opts.frameRate;
        open(writer);
        fprintf('trackMovie: writing to %s ...\n', opts.outputFile);
    end

    %% Main animation loop.
    trailLen = opts.trailLength;

    for iT = 1:nT

        %% Background frame.
        bg = imF(:, :, iT);
        rgbFrame = repmat(bg, [1, 1, 3]);

        %% Mito outlines (white).
        ol = outlines{iT};
        if ~isempty(ol)
            for c = 1:3
                ch = rgbFrame(:,:,c);
                ch(ol) = 1;
                rgbFrame(:,:,c) = ch;
            end
        end

        %% Draw track trails and current positions.
        rgbFrame = drawTracks(rgbFrame, tracks, trackFrameMap, iT, ...
            trailLen, trackColors, opts, speedMax, p.streamingThreshold);

        %% Display.
        if iT == 1
            hImg = imshow(rgbFrame, 'Parent', ax);
            hold(ax, 'on');
            % Frame counter text.
            hText = text(ax, 5, 12, frameLabel(iT, nT, p), ...
                'Color', 'w', 'FontSize', 9, 'FontName', 'monospaced', ...
                'VerticalAlignment', 'top');
        else
            set(hImg, 'CData', rgbFrame);
            set(hText, 'String', frameLabel(iT, nT, p));
        end

        %% Track ID labels (drawn fresh each frame to handle movement).
        clearIDLabels(ax);
        if opts.showIDs
            drawIDLabels(ax, tracks, trackFrameMap, iT, trackColors);
        end

        drawnow;

        if ~isempty(writer)
            writeVideo(writer, getframe(fig));
        else
            pause(opts.framePause);
        end
    end

    if ~isempty(writer)
        close(writer);
        fprintf('trackMovie: saved %d frames to %s\n', nT, opts.outputFile);
    end
end

% =========================================================================
% Local helpers
% =========================================================================

function opts = parseOptions(args)
    opts.trailLength  = 15;
    opts.colorBy      = 'track';
    opts.showIDs      = true;
    opts.frameRate    = 10;
    opts.outputFile   = '';
    opts.framePause   = 0.05;

    i = 1;
    while i <= numel(args) - 1
        key = args{i};
        val = args{i+1};
        switch lower(key)
            case 'traillength',  opts.trailLength  = val;
            case 'colorby',      opts.colorBy      = val;
            case 'showids',      opts.showIDs      = val;
            case 'framerate',    opts.frameRate     = val;
            case 'outputfile',   opts.outputFile    = val;
            case 'framepause',   opts.framePause    = val;
            otherwise
                warning('trackMovie: unknown option ''%s'' ignored.', key);
        end
        i = i + 2;
    end
end

% -------------------------------------------------------------------------

function outlines = precomputeOutlines(labelStack, nT)
    outlines = cell(nT, 1);
    for iT = 1:nT
        bw = labelStack(:,:,iT) > 0;
        if any(bw(:))
            outlines{iT} = bw & ~imerode(bw, ones(3,3));
        end
    end
end

% -------------------------------------------------------------------------

function map = buildFrameMap(tracks, nT)
    % map{iT} = struct array with one entry per track visible at frame iT.
    % Each entry: .trackIdx, .x, .y, .speedUm, .streaming
    map = cell(nT, 1);
    for i = 1:numel(tracks)
        tr = tracks(i);
        for k = 1:numel(tr.frames)
            f = tr.frames(k);
            if f < 1 || f > nT, continue, end
            entry.trackIdx  = i;
            entry.x         = tr.x(k);
            entry.y         = tr.y(k);
            entry.speedUm   = tr.speedUm(k);
            entry.streaming = tr.streaming(k);
            map{f} = [map{f}, entry]; %#ok<AGROW>
        end
    end
end

% -------------------------------------------------------------------------

function rgb = drawTracks(rgb, tracks, ~, iT, trailLen, colors, opts, speedMax, streamThresh)
    % Draw filled trail lines for each track visible up to iT.
    nTracks = numel(tracks);
    for i = 1:nTracks
        tr   = tracks(i);
        col  = colors(i, :);

        % Find observations within the trail window ending at iT.
        inWindow = tr.frames >= (iT - trailLen) & tr.frames <= iT;
        if ~any(inWindow), continue, end

        fx = tr.x(inWindow);
        fy = tr.y(inWindow);
        fs = tr.streaming(inWindow);
        fsp= tr.speedUm(inWindow);

        if numel(fx) < 2, continue, end

        % Draw segment-by-segment so streaming can be coloured red.
        for s = 1:numel(fx)-1
            x1 = fx(s);  y1 = fy(s);
            x2 = fx(s+1); y2 = fy(s+1);

            if strcmp(opts.colorBy, 'speed')
                % Interpolate speed to a blue-red colour.
                spVal = max(fsp(s), 0) / speedMax;
                segCol = speedColor(spVal);
            else
                segCol = col;
            end

            % Streaming override: red, thicker (approximated by 3-px line).
            if fs(s) || fs(s+1)
                segCol = [1, 0.15, 0.15];
                rgb = drawLine(rgb, x1, y1, x2, y2, segCol, 2);
            else
                rgb = drawLine(rgb, x1, y1, x2, y2, segCol, 1);
            end
        end

        % Current position: filled disc.
        lastIdx = find(tr.frames == iT, 1);
        if ~isempty(lastIdx)
            cx = round(tr.x(lastIdx));
            cy = round(tr.y(lastIdx));
            rgb = drawDisc(rgb, cx, cy, 3, col);
        end
    end
end

% -------------------------------------------------------------------------

function rgb = drawLine(rgb, x1, y1, x2, y2, col, thickness)
    % Bresenham-style anti-aliased line drawing into an RGB image.
    [nY, nX, ~] = size(rgb);
    nPts = max(abs(x2-x1), abs(y2-y1)) * 3 + 2;
    xs = round(linspace(x1, x2, nPts));
    ys = round(linspace(y1, y2, nPts));

    % Thicken by painting a small neighbourhood.
    offsets = 0:thickness-1;
    offsets = offsets - floor(thickness/2);

    for k = 1:numel(xs)
        for od = offsets
            px = xs(k) + od;
            py = ys(k) + od;
            if px < 1 || px > nX || py < 1 || py > nY, continue, end
            rgb(py, px, 1) = col(1);
            rgb(py, px, 2) = col(2);
            rgb(py, px, 3) = col(3);
            % Also paint orthogonal offset for thicker lines.
            px2 = xs(k) - od;
            py2 = ys(k) + od;
            if px2 >= 1 && px2 <= nX && py2 >= 1 && py2 <= nY
                rgb(py2, px2, 1) = col(1);
                rgb(py2, px2, 2) = col(2);
                rgb(py2, px2, 3) = col(3);
            end
        end
    end
end

% -------------------------------------------------------------------------

function rgb = drawDisc(rgb, cx, cy, r, col)
    [nY, nX, ~] = size(rgb);
    for dy = -r:r
        for dx = -r:r
            if dx^2 + dy^2 <= r^2
                px = cx + dx;
                py = cy + dy;
                if px >= 1 && px <= nX && py >= 1 && py <= nY
                    rgb(py, px, 1) = col(1);
                    rgb(py, px, 2) = col(2);
                    rgb(py, px, 3) = col(3);
                end
            end
        end
    end
end

% -------------------------------------------------------------------------

function col = speedColor(frac)
    % Map frac in [0,1] to blue(0) -> cyan -> yellow -> red(1).
    frac = max(0, min(1, frac));
    cmap = [0.0, 0.2, 1.0;   % blue
            0.0, 0.8, 0.8;   % cyan
            1.0, 0.9, 0.0;   % yellow
            1.0, 0.1, 0.1];  % red
    t   = frac * (size(cmap,1) - 1) + 1;
    lo  = floor(t);  hi = min(lo + 1, size(cmap,1));
    w   = t - lo;
    col = cmap(lo,:) * (1-w) + cmap(hi,:) * w;
end

% -------------------------------------------------------------------------

function clearIDLabels(ax)
    % Remove text objects tagged as track ID labels.
    kids = findobj(ax, 'Type', 'text', 'Tag', 'trackID');
    delete(kids);
end

function drawIDLabels(ax, tracks, trackFrameMap, iT, colors)
    visible = trackFrameMap{iT};
    for k = 1:numel(visible)
        i   = visible(k).trackIdx;
        col = colors(i, :);
        text(ax, visible(k).x + 4, visible(k).y - 4, ...
             num2str(tracks(i).id), ...
             'Color', col, 'FontSize', 6, ...
             'VerticalAlignment', 'bottom', ...
             'Tag', 'trackID');
    end
end

% -------------------------------------------------------------------------

function s = frameLabel(iT, nT, p)
    tSec = (iT - 1) * p.frameInterval;
    s = sprintf('Frame %3d / %d   t = %.1f s', iT, nT, tSec);
end
