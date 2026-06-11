%DEMOMITOSIMTRACKING  End-to-end demo: generate synthetic movie, track, evaluate.
%
%   Run this script from the OrganelleTracking_sandbox root (or any folder
%   on the MATLAB path that includes src/, utils/, demos/).
%
%   Steps
%   -----
%   1. Generate a 512x512 x 50-frame synthetic movie with 'mixed' motion.
%   2. Run trackCompute on a binarised / labelled version of the movie.
%   3. Play the movie with ground-truth overlays.
%   4. Show the demoTrackMethods summary figure.
%   5. Evaluate tracking accuracy against the simulator ground truth.
%   6. Save the raw frames as a TIFF stack.
%
%   The script assumes the Computer Vision Toolbox is available
%   (assignDetectionsToTracks) and Image Processing Toolbox (regionprops,
%   imfilter, fspecial, poissrnd).

%% 0. Parameters -------------------------------------------------------

simP = mitoSimParamsDefault();
simP.numMitochondria = 20;
simP.numFrames       = 50;
simP.imageSize       = [512, 512];
simP.motionType      = 'mixed';
simP.poissonNoise    = true;
simP.gaussianNoise   = 15;

%% 1. Generate --------------------------------------------------------

[frames, groundTruth] = mitoSimGenerate(simP);

%% 2. Segment (simple threshold) -> label stack -----------------------
% In practice you would feed the output of the AnalyzER organelle
% segmentation pipeline here.  For this demo a fixed threshold suffices.

fprintf('Segmenting frames...\n');
threshold  = simP.background * 1.5;
labelStack = zeros(size(frames), 'uint16');

for t = 1:size(frames, 3)
    bw              = frames(:,:,t) > threshold;
    bw              = imopen(bw, strel('disk', 1));
    labelStack(:,:,t) = uint16(bwlabel(bw));
end

%% 3. Track -----------------------------------------------------------

trkP                 = trackParamsDefault();
trkP.pixelSize       = simP.pixelSize;
trkP.frameInterval   = simP.frameInterval;
trkP.maxDisplacement = 20;
trkP.minTrackLength  = 5;

tracks = trackCompute(labelStack, trkP);
fprintf('Tracking found %d tracks.\n', numel(tracks));

%% 4. Play movie (raw + GT overlay) -----------------------------------

mitoSimPlayMovie(frames, groundTruth, true);

%% 5. Summary figure + animated movie ---------------------------------

% Use frame 1 as the raw image for demoTrackMethods.
imRaw = double(frames(:,:,1));

demoTrackMethods(labelStack, imRaw, trkP);

trackMovie(labelStack, imRaw, tracks, trkP, ...
    'TrailLength', 10, 'ColorBy', 'track', 'FramePause', 0.04);

%% 6. Evaluate --------------------------------------------------------

metrics = mitoSimEvaluateTracking(tracks, groundTruth, trkP.minTrackLength);

%% 7. Save TIFF stack -------------------------------------------------

mitoSimSaveMovie('synthetic_mitochondria.tif', frames);
