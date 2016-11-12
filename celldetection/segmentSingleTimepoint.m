function [stats,params,bw] = segmentSingleTimepoint(I,bw,params)

scndtry = 0;
status = 0;

bw2 = bw;
while scndtry <= 2
    if status == 0
        [stats,status] = measureStats(I,bw2,params,0);
    elseif status == 2
        [~,bw2] = watershedCellProfiler(I,bw,'distance','hough','80,1000',0,0,0);
        [stats,status] = measureStats(I,bw2,params,0);
        disp('do watershedding...')
        scndtry = scndtry +1;
    elseif status == 3
        bw2 = imopen(bw,strel('disk',2));
        [stats,status] = measureStats(I,bw2,params,0);
        scndtry = scndtry +1;
    elseif status == 4
        scndtry = 3;
        bw2 = bw;
        [stats,status] = measureStats(I,bw2,params,1);
        disp('use bw...')
    elseif status == 1
        scndtry = 3;
    else
        error ('no status')
    end
    
    if numel(fields(stats)) == 0
        stats = struct();
    end
    
end


params.prevCell = params.currentCell;
if numel(fields(stats)) == 0
        params.lastCentroid = params.center;
        params.prevSize = nan;
elseif ~isnan(stats.Area)
        params.lastCentroid = params.center;
        params.prevSize = nan;
else
    disp(stats.Centroid)
    disp(stats.Area)
    toc
    params.prevSize = stats.Area;
    bw = bw2;
end
end

