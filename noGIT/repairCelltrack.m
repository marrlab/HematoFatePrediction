load('F:\Data\Philipp\FatePrediction_2013\NatMethDemo\demoTracks.mat')

unicells = unique(track.cellNr);
trackAttributes = fields(track);
moviepath = 'E:\140206PH8\';
moviename = '140206PH8';
for cellnr = unicells
    celltrack = struct();
    
    for f =1:numel(trackAttributes)
        celltrack.(trackAttributes{f}) = track.(trackAttributes{f})(track.cellNr == cellnr);
    end
    
    
    celltrack_corrected = struct();
    
    alltps = min(celltrack.timepoint):max(celltrack.timepoint);
    [~,sortedtpidx] = sort(celltrack.timepoint);
    
    celltrack_corrected.missingtps = ~ismember(alltps,celltrack.timepoint);
    
    if ~isempty(celltrack_corrected.missingtps)
        celltrackpos = unique(celltrack.positionIndex);
        moviestart = getmoviestart(moviepath);
        poslog = cell(numel(celltrackpos),1);
        for p=celltrackpos
            poslog{p==celltrackpos} = positionLogFileReader(...
                sprintf('%s/%s_p%.4i/%s_p%.4i.log',moviepath,moviename,p,moviename,p));
        end
    end
        
    goodidx = find(ismember(alltps,celltrack.timepoint));
    badidx = find(~ismember(alltps,celltrack.timepoint));
        
    fieldz = fields(celltrack);
    for f=fieldz'
        if strcmp(f{:},'timepoint') || strcmp(f{:},'positionIndex')  || strcmp(f{:},'absoluteTime')
            celltrack_corrected.timepoint = alltps;
            celltrack_corrected.absoluteTime = nan(1,numel(alltps));
            celltrack_corrected.absoluteTime(goodidx) = celltrack.absoluteTime(sortedtpidx);
            if ~isempty(celltrack_corrected.missingtps)
                celltrack_corrected.positionIndex = nan(1,numel(alltps));
                celltrack_corrected.positionIndex(goodidx) = celltrack.positionIndex(sortedtpidx);
                interpvals = interp1(goodidx,celltrack_corrected.positionIndex(goodidx),badidx,'nearests','extrap');
                celltrack_corrected.positionIndex(badidx) = interpvals;
                
                for i = badidx
                    logtpidx = poslog{celltrackpos == celltrack_corrected.positionIndex(i)}.timepoint == alltps(i) &...
                        poslog{celltrackpos == celltrack_corrected.positionIndex(i)}.wavelength == bfwl;
                    celltrack_corrected.absoluteTime(i) = poslog{celltrackpos == celltrack_corrected.positionIndex(i)}.absoluteTime(logtpidx)-moviestart;
                end
                
            end
        elseif strcmp(f{:},'filename') || strcmp(f{:},'annotator')
            celltrack_corrected.(f{:}) = cell(1,numel(alltps));
            celltrack_corrected.(f{:})(goodidx) = celltrack.(f{:})(sortedtpidx);
        else
            celltrack_corrected.(f{:}) = nan(1,numel(alltps));
            celltrack_corrected.(f{:})(goodidx) = celltrack.(f{:})(sortedtpidx);
            interpvals = interp1(goodidx,celltrack_corrected.(f{:})(goodidx),badidx,'nearests','extrap');
            celltrack_corrected.(f{:})(badidx) = interpvals;
        end
    end
    correctedTrackAttributes = fields(celltrack_corrected);
    if find(cellnr == unicells) == 1
        tracks_corrected = celltrack_corrected;
    else
        for f =1:numel(correctedTrackAttributes)
            tracks_corrected.(correctedTrackAttributes{f}) = [tracks_corrected.(correctedTrackAttributes{f}), celltrack_corrected.(correctedTrackAttributes{f})];
        end
    end
end
track = tracks_corrected;
save('F:\Data\Philipp\FatePrediction_2013\NatMethDemo\demoTracks_corrected.mat','track')