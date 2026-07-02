function detections = trackDetectionsExtract(labelStack, p, wb, wbLo, wbHi)
%TRACKDETECTIONSEXTRACT  Extract per-frame shape features from a label stack.
%
%   detections = trackDetectionsExtract(labelStack, p)
%   detections = trackDetectionsExtract(labelStack, p, wb, wbLo, wbHi)
%
%   labelStack : [nY nX nT] uint16 label matrix.  Each non-zero value is an
%                object label for that frame.  Labels need not be consistent
%                across frames -- tracking handles that.
%   p          : parameter struct. Optional field p.boundaryLabels: a
%                {nT x 1} (or 1 x nT) cell array, one entry per frame, each
%                holding the vector of label IDs touching the cell
%                boundary in that frame (see organelleBoundaryFlag /
%                app.variables.organelle.refine.boundaryLabels). Missing or
%                empty -> touchesBoundary is false for every detection.
%   wb, wbLo, wbHi : optional waitbar handle and fractional range [wbLo,wbHi]
%                    to update within.  Caller owns the waitbar lifecycle.
%                    Cancel is signalled by setappdata(wb,'cancel',true).
%
%   detections : [nT x 1] cell array.  Each cell contains a struct with
%                fields (all column vectors, one row per object):
%
%       x, y          Centroid (px, 1-based, column-major MATLAB convention).
%       area          Object area (px^2).
%       majorAxis     Major axis length (px).
%       minorAxis     Minor axis length (px).
%       aspectRatio   majorAxis / minorAxis  (>=1; 1 = circle).
%       orientation   Axis angle (deg, range -90..90, regionprops convention).
%       eccentricity  0 = circle, 1 = line segment.
%       labelID       Label value in the input stack for this object.
%       touchesBoundary  Logical; true if this object touches the cell
%                        boundary in this frame (see p.boundaryLabels).
%
%   Empty cell entries (no objects in that frame) contain [].

    useWb = nargin >= 3 && ~isempty(wb) && ishandle(wb);

    boundaryLabels = {};
    if nargin >= 2 && isstruct(p) && isfield(p,'boundaryLabels') && ~isempty(p.boundaryLabels)
        boundaryLabels = p.boundaryLabels;
    end

    [~, ~, nT] = size(labelStack);
    detections = cell(nT, 1);

    propNames = {'Centroid','Area','MajorAxisLength','MinorAxisLength', ...
                 'Orientation','Eccentricity'};

    for iT = 1:nT

        if useWb
            frac = wbLo + (iT / nT) * (wbHi - wbLo);
            waitbar(frac, wb, sprintf('Extracting features: frame %d / %d', iT, nT));
            if getappdata(wb, 'cancel')
                error('trackCompute:cancelled', ...
                    'Tracking cancelled by user at feature extraction frame %d / %d.', iT, nT);
            end
        end

        frame = labelStack(:, :, iT);
        if ~any(frame(:))
            continue
        end
        S = regionprops(frame, propNames{:});
        if isempty(S)
            continue
        end

        xy  = vertcat(S.Centroid);          % n x 2
        maj = [S.MajorAxisLength]';
        mn  = [S.MinorAxisLength]';

        d.x            = xy(:, 1);
        d.y            = xy(:, 2);
        d.area         = [S.Area]';
        d.majorAxis    = maj;
        d.minorAxis    = mn;
        d.aspectRatio  = maj ./ max(mn, 1);  % guard zero minor axis
        d.orientation  = [S.Orientation]';
        d.eccentricity = [S.Eccentricity]';
        d.labelID      = (1 : numel(S))';

        d.touchesBoundary = false(numel(S), 1);
        if iT <= numel(boundaryLabels) && ~isempty(boundaryLabels{iT})
            d.touchesBoundary = ismember(d.labelID, boundaryLabels{iT});
        end

        detections{iT} = d;
    end
end
