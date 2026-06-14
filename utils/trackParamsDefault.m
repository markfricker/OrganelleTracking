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
%   costNonAssignment  Cost of starting a new track or ending an existing
%                      one.  Set higher to prefer gaps over wrong links.
%
%   Cost weights  (weighted sum, each term in [0,1]^2)
%   wDisp    Normalised squared displacement to predicted position.
%   wOrient  Normalised squared axis-angle difference (mod 180 deg).
%   wShape   Normalised squared aspect-ratio difference.
%   wArea    Normalised squared relative area change (0 = off).
%   wDir     Squared motion-direction inconsistency (0 = off).
%
%   Linking gates (toggles)
%   -----------------------
%   usePredictedGate   true  Gate candidate search on predicted position
%                            rather than last known position.  Tighter
%                            and more appropriate for fast-moving objects;
%                            the gate radius is maxDisplacement from
%                            predX/predY.
%                      false Legacy behaviour: gate is 2*maxDisplacement
%                            from last known position (wider, less strict).
%
%   useAreaCost        true  Penalise large area changes between frames.
%                            Helps reject incorrect links when a slow mito
%                            is nearby and has a different size.
%
%   useDirCost         true  Penalise candidates whose displacement from
%                            the predicted position is inconsistent with
%                            the track velocity direction.  Most effective
%                            for streaming mitochondria with established
%                            velocity estimates.
%
%   Velocity prediction
%   -------------------
%   velocityAlpha      EMA weight for velocity update (0=instant, 1=frozen).
%   velocityMinFrames  Frames before velocity prediction is trusted.
%
%   Gap closing
%   -----------
%   gapMax         Maximum missing frames that can be bridged.
%   gapCostScale   Cost multiplier per missing frame.
%
%   Post-filter
%   -----------
%   minTrackLength  Tracks shorter than this (frames) are discarded.
%
%   Streaming detection
%   -------------------
%   streamingThreshold  Speed (px/frame) above which a frame is flagged
%                       as streaming.
%
%   Acquisition
%   -----------
%   pixelSize      um per pixel.
%   frameInterval  Seconds per frame.

p.maxDisplacement   = 20;     % px  (~1.8 um at 90 nm/px)
p.costNonAssignment = 50;

p.wDisp             = 1.0;
p.wOrient           = 0.3;
p.wShape            = 0.1;
p.wArea             = 0.2;    % set to 0 to disable area cost
p.wDir              = 0.3;    % set to 0 to disable direction cost

p.usePredictedGate  = true;   % gate on predicted position (recommended)
p.useAreaCost       = true;
p.useDirCost        = true;
p.useFwdBwd         = false;  % run backward pass and tag unconfirmed links

p.velocityAlpha     = 0.7;
p.velocityMinFrames = 2;

p.gapMax            = 3;
p.gapCostScale      = 1.5;

p.minTrackLength    = 5;

p.streamingThreshold = 11;    % px/frame (~1 um/s at 90 nm/px, 0.4 s/frame)

p.pixelSize         = 0.09;   % um/px
p.frameInterval     = 0.4;    % s/frame
end
