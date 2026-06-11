function [frames, groundTruth] = mitoSimGenerate(p)
%MITOSIMGENERATE  Generate a synthetic mitochondria time series with ground truth.
%
%   [frames, groundTruth] = mitoSimGenerate(p)
%
%   Produces a confocal-like image stack with ovoid mitochondria undergoing
%   configurable motion (Brownian, directed, or mixed), PSF blurring, Poisson
%   noise, and slow photobleaching.
%
%   INPUT
%     p : parameter struct from mitoSimParamsDefault (missing fields are filled).
%
%   OUTPUTS
%     frames      : [nY nX nT] uint16 image stack.
%     groundTruth : [nT x 1] cell array.  Each cell contains a struct array
%                   with one element per mitochondrion:
%                     id, x, y, vx, vy, majorAxis, minorAxis, theta, intensity
%
%   EXAMPLE
%     p = mitoSimParamsDefault();
%     p.numMitochondria = 20;
%     p.numFrames       = 50;
%     p.motionType      = 'mixed';
%     [frames, gt] = mitoSimGenerate(p);
%     mitoSimPlayMovie(frames, gt, p);
%
%   NOTE
%     Requires Image Processing Toolbox (imfilter, fspecial, poissrnd).

    p = fillDefaults(p);

    fprintf('Generating synthetic mitochondria movie...\n');
    fprintf('  Image: %d x %d px  |  %d mitochondria  |  %d frames  |  motion: %s\n', ...
        p.imageSize(2), p.imageSize(1), p.numMitochondria, p.numFrames, p.motionType);

    % --- PSF ------------------------------------------------------------------
    psf = generatePSF(p);

    % --- Initialise objects ---------------------------------------------------
    mito = initMitochondria(p);
    fprintf('  Placed %d / %d mitochondria.\n', numel(mito), p.numMitochondria);

    % --- Render frames --------------------------------------------------------
    frames      = zeros([p.imageSize, p.numFrames], 'uint16');
    groundTruth = cell(p.numFrames, 1);

    fprintf('  Rendering: ');
    for t = 1:p.numFrames
        if mod(t, 10) == 0
            fprintf('%d ', t);
        end
        frames(:,:,t) = uint16(renderFrame(mito, psf, p));
        groundTruth{t} = mito;
        if t < p.numFrames
            mito = updatePositions(mito, p);   % returns updated struct (fixes pass-by-value)
        end
    end
    fprintf('done.\n');
end

% =========================================================================
% Local: PSF
% =========================================================================

function psf = generatePSF(p)
    % Airy-disk approximation: PSF radius = 0.61 * lambda / NA
    r_airy    = 0.61 * p.psfLambda / p.psfNA;           % um
    sigma_psf = r_airy / (2.355 * p.pixelSize);          % px (FWHM -> sigma)

    psfSize = ceil(6 * sigma_psf);
    if mod(psfSize, 2) == 0
        psfSize = psfSize + 1;
    end

    psf = fspecial('gaussian', psfSize, sigma_psf);
    psf = psf / sum(psf(:));

    fprintf('  PSF: sigma = %.2f px (%.3f um)\n', sigma_psf, sigma_psf * p.pixelSize);
end

% =========================================================================
% Local: initialise non-overlapping mitochondria
% =========================================================================

function mito = initMitochondria(p)
    mito        = [];
    maxAttempts = 1000;
    margin      = 10;   % px from edge

    for i = 1:p.numMitochondria
        placed   = false;
        attempts = 0;

        while ~placed && attempts < maxAttempts
            attempts = attempts + 1;

            majorAxis = p.majorAxisRange(1) + diff(p.majorAxisRange) * rand();
            minorAxis = p.minorAxisRange(1) + diff(p.minorAxisRange) * rand();
            maxRadius = max(majorAxis, minorAxis);

            x = margin + maxRadius + ...
                (p.imageSize(2) - 2*margin - 2*maxRadius) * rand();
            y = margin + maxRadius + ...
                (p.imageSize(1) - 2*margin - 2*maxRadius) * rand();

            theta     = 2 * pi * rand();
            intensity = p.intensityRange(1) + diff(p.intensityRange) * rand();

            [vx, vy] = sampleVelocity(p);

            % Check overlap with already-placed objects.
            overlap = false;
            for j = 1:numel(mito)
                dist    = hypot(x - mito(j).x, y - mito(j).y);
                minDist = max(majorAxis, mito(j).majorAxis) + ...
                          max(minorAxis, mito(j).minorAxis);
                if dist < minDist
                    overlap = true;
                    break
                end
            end

            if ~overlap
                placed = true;
            end
        end

        if placed
            ns.id        = i;
            ns.x         = x;
            ns.y         = y;
            ns.vx        = vx;
            ns.vy        = vy;
            ns.majorAxis = majorAxis;
            ns.minorAxis = minorAxis;
            ns.theta     = theta;
            ns.intensity = intensity;
            if isempty(mito)
                mito = ns;
            else
                mito(end+1) = ns; %#ok<AGROW>
            end
        else
            warning('mitoSimGenerate:placementFailed', ...
                'Could not place mitochondrion %d after %d attempts.', i, maxAttempts);
        end
    end
end

function [vx, vy] = sampleVelocity(p)
    switch p.motionType
        case 'directed'
            angle = 2 * pi * rand();
            speed = p.maxVelocity * rand() / p.pixelSize;  % um/s -> px/s
            vx    = speed * cos(angle);
            vy    = speed * sin(angle);
        case 'mixed'
            if rand() > 0.5
                vx = randn() * 0.5;
                vy = randn() * 0.5;
            else
                angle = 2 * pi * rand();
                speed = p.maxVelocity * rand() / p.pixelSize;
                vx    = speed * cos(angle);
                vy    = speed * sin(angle);
            end
        otherwise  % 'brownian'
            vx = randn() * 0.5;
            vy = randn() * 0.5;
    end
end

% =========================================================================
% Local: render one frame
% =========================================================================

function frame = renderFrame(mito, psf, p)
    frame = zeros(p.imageSize);

    for i = 1:numel(mito)
        frame = frame + renderMitochondrion(mito(i), p);
    end

    % PSF convolution.
    frame = imfilter(frame, psf, 'replicate');

    % Background offset.
    frame = frame + p.background;

    % Poisson noise: scale to photon counts, sample, scale back.
    if p.poissonNoise
        above = frame(frame > p.background);
        if ~isempty(above)
            scale = p.photonCount / mean(above);
            frame = double(poissrnd(frame * scale)) / scale;
        end
    end

    % Additive Gaussian noise.
    if p.gaussianNoise > 0
        frame = frame + p.gaussianNoise * randn(p.imageSize);
    end

    frame(frame < 0) = 0;
end

function img = renderMitochondrion(m, p)
    % Elliptical Gaussian blob at position (m.x, m.y).
    [X, Y]   = meshgrid(1:p.imageSize(2), 1:p.imageSize(1));
    Xr = (X - m.x) .* cos(m.theta) + (Y - m.y) .* sin(m.theta);
    Yr = -(X - m.x) .* sin(m.theta) + (Y - m.y) .* cos(m.theta);

    sigX = m.majorAxis / 2;
    sigY = m.minorAxis / 2;

    img = m.intensity .* exp(-((Xr / sigX).^2 + (Yr / sigY).^2) / 2);
end

% =========================================================================
% Local: advance positions one frame
% =========================================================================

function mito = updatePositions(mito, p)
    % Returns updated struct array.  Must return to avoid pass-by-value loss.
    dt     = p.frameInterval;
    margin = 10;

    for i = 1:numel(mito)
        % Brownian perturbation (always present).
        bx = randn() * sqrt(dt) * 0.3;
        by = randn() * sqrt(dt) * 0.3;

        mito(i).x = mito(i).x + mito(i).vx * dt + bx;
        mito(i).y = mito(i).y + mito(i).vy * dt + by;

        % Elastic wall bounce with damping.
        maxR = max(mito(i).majorAxis, mito(i).minorAxis);

        if mito(i).x < margin + maxR
            mito(i).x  = margin + maxR;
            mito(i).vx = -mito(i).vx * 0.8;
        elseif mito(i).x > p.imageSize(2) - margin - maxR
            mito(i).x  = p.imageSize(2) - margin - maxR;
            mito(i).vx = -mito(i).vx * 0.8;
        end

        if mito(i).y < margin + maxR
            mito(i).y  = margin + maxR;
            mito(i).vy = -mito(i).vy * 0.8;
        elseif mito(i).y > p.imageSize(1) - margin - maxR
            mito(i).y  = p.imageSize(1) - margin - maxR;
            mito(i).vy = -mito(i).vy * 0.8;
        end

        % Slow tumble.
        mito(i).theta = mito(i).theta + 0.05 * randn();

        % Slow photobleaching + intensity jitter.
        mito(i).intensity = mito(i).intensity * 0.998 + randn() * 10;
    end
end

% =========================================================================
% Local: fill missing parameter fields from defaults
% =========================================================================

function p = fillDefaults(p)
    d   = mitoSimParamsDefault();
    fns = fieldnames(d);
    for k = 1:numel(fns)
        if ~isfield(p, fns{k})
            p.(fns{k}) = d.(fns{k});
        end
    end
end
