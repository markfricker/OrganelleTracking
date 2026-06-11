function mitoSimSaveMovie(filename, frames)
%MITOSIMSAMOVIE  Save a synthetic mitochondria stack as a multi-page TIFF.
%
%   mitoSimSaveMovie(filename, frames)
%
%   INPUTS
%     filename : output file path (e.g. 'synthetic_mito.tif').
%     frames   : [nY nX nT] uint16 image stack (from mitoSimGenerate).
%
%   The file is written uncompressed so it can be opened directly in FIJI /
%   ImageJ as a hyperstack.

    nT = size(frames, 3);
    fprintf('Saving %d frames to %s ... ', nT, filename);

    for t = 1:nT
        if t == 1
            imwrite(frames(:,:,t), filename, 'tiff', ...
                'Compression', 'none', 'WriteMode', 'overwrite');
        else
            imwrite(frames(:,:,t), filename, 'tiff', ...
                'Compression', 'none', 'WriteMode', 'append');
        end
    end

    fprintf('done.\n');
end
