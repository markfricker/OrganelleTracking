function p = trackParamsDefault()
%TRACKPARAMSDEFAULT  Factory default parameters for organelle tracking.
%
%   p = trackParamsDefault()
%
%   All spatial parameters are in pixels; convert using p.pixelSize (um/px)
%   and p.frameInterval (s) when computing physical quantities.
%
%   Linking parameters
%   ------------------
%   maxDisplacement    Maximum allowed frame-to-frame displacement (px).
%                      Detections further away than this are never linked
%                      regardless of cost.
%   costNonAssignment  Cost of starting a new track or ending an existing
%                      one.  Set to (maxDisplacement/2)^2 so a displacement
%                      beyond half the search radius prefers birth/death.
%
%   Cost weights  (wDisp + wOrient + wShape = cost range [0, sum])
%   wDisp          Weight on normalised squared displacement.
%   wOrient        Weight on normalised squared orientation difference.
%   wShape         Weight on normalised squared aspect-ratio difference.
%
%   Velocity prediction
%   -------------------
%   velocityAlpha  Exponential moving average weight for velocity update.
%                  0 = always use instantaneous velocity (noisy),
%                  1 = never update (fixed from first frame).
%                  0.7 gives a ~3-frame smoothing window.
%   velocityMinFrames  Minimum track length before velocity prediction is
%                      applied.  Short tracks use last-position linking.
%
%   Gap closing
%   -----------
%   gapMax         Maximum number of missing frames that can be bridged.
%   gapCostScale   Cost multiplier per missing frame (discourages long gaps).
%
%   Post-filter
%   -----------
%   minTrackLength  Tracks shorter than this (frames) are discarded.
%
%   Acquisition
%   -----------
%   pixelSize      um per pixel (for physical-unit outputs).
%   frameInterval  Seconds per frame (for velocity in um/s).

p.maxDisplacement   = 20;     % px  (~1.8 um at 90 nm/px)
p.costNonAssignment = 50;     % dimensionless (see above)

p.wDisp             = 1.0;
p.wOrient           = 0.3;
p.wShape            = 0.1;

p.velocityAlpha     = 0.7;    % EMA smoothing (0 = instant, 1 = frozen)
p.velocityMinFrames = 2;      % frames before prediction is trusted

p.gapMax            = 3;      % frames
p.gapCostScale      = 1.5;    % cost x scale per gap frame

p.minTrackLength    = 5;      % frames

p.pixelSize         = 0.09;   % um/px
p.frameInterval     = 0.4;    % s/frame
end
