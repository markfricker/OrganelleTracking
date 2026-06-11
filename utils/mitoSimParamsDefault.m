function p = mitoSimParamsDefault()
%MITOSIMPARAMSDEFAULT  Factory parameters for the synthetic mitochondria simulator.
%
%   p = mitoSimParamsDefault()
%
%   Returns a parameter struct with sensible defaults matching a 512x512
%   confocal acquisition at 0.1 um/px and 0.4 s/frame (GFP channel).
%
%   FIELDS
%   Image
%     imageSize         [nY nX] in pixels            default [512 512]
%     pixelSize         um/pixel                     default 0.1
%
%   Mitochondria
%     numMitochondria   number of objects             default 25
%     majorAxisRange    [min max] px                  default [3 8]
%     minorAxisRange    [min max] px                  default [2 4]
%     intensityRange    [min max] mean intensity      default [800 1200]
%
%   Motion
%     maxVelocity       um/s  (directed/mixed mode)   default 4
%     frameInterval     s/frame                       default 0.4
%     numFrames         total frames to simulate      default 100
%     motionType        'brownian' | 'directed' | 'mixed'   default 'brownian'
%
%   Imaging / noise
%     psfNA             objective numerical aperture  default 1.4
%     psfLambda         emission wavelength um (GFP)  default 0.520
%     background        background intensity offset   default 100
%     poissonNoise      logical -- add Poisson noise  default true
%     gaussianNoise     std of additive Gaussian      default 20
%     photonCount       avg photons / mitochondrion   default 1000

    % Image
    p.imageSize         = [512, 512];
    p.pixelSize         = 0.1;          % um/px

    % Mitochondria
    p.numMitochondria   = 25;
    p.majorAxisRange    = [3, 8];       % px (0.3 - 0.8 um at 0.1 um/px)
    p.minorAxisRange    = [2, 4];       % px (0.2 - 0.4 um)
    p.intensityRange    = [800, 1200];  % mean intensity counts

    % Motion
    p.maxVelocity       = 4;            % um/s
    p.frameInterval     = 0.4;          % s/frame
    p.numFrames         = 100;
    p.motionType        = 'brownian';   % 'brownian', 'directed', 'mixed'

    % Imaging / noise
    p.psfNA             = 1.4;
    p.psfLambda         = 0.520;        % um (GFP emission)
    p.background        = 100;
    p.poissonNoise      = true;
    p.gaussianNoise     = 20;           % std of additive Gaussian noise
    p.photonCount       = 1000;         % avg photons per mitochondrion
end
