function [tgtLabel, found] = trackOriginalLink(baseTracks, frame, srcLabel)
%TRACKORIGINALLINK  What label did the tracker originally link srcLabel@frame to?
%
%   [tgtLabel, found] = trackOriginalLink(baseTracks, frame, srcLabel)
%
% Looks up the track that owns srcLabel at `frame` in the ORIGINAL
% (unedited) baseTracks, and returns whatever label it continued to at
% frame+1, if any. found=false if srcLabel has no track at `frame`, the
% track ends there (no frame+1 observation), or frame+1 is a gap-
% interpolated placeholder (labelID==0).

tgtLabel = 0;  found = false;
for k = 1:numel(baseTracks)
    fIdx = find(baseTracks(k).frames == frame, 1);
    if isempty(fIdx) || baseTracks(k).labelID(fIdx) ~= srcLabel, continue; end
    fIdx2 = find(baseTracks(k).frames == frame+1, 1);
    if ~isempty(fIdx2) && baseTracks(k).labelID(fIdx2) > 0
        tgtLabel = baseTracks(k).labelID(fIdx2);
        found = true;
    end
    return
end
end
