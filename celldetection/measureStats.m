%tries to measure the segmented objects. Object status is one of the
%following status numbers:
% 0 = finally nothing found
% 1 = found a cell which has a reasonable size and shape
% 2 = found a cell with good size but problematic shape --> apply
%     watershedding
% 3 = found object which is too small --> try to mask beads and retry
%     segmentation
% 4 = nothing found

function [cellStats,status] = measureStatsPH(I_c,L,params,usebw)
stats=regionprops(L,I_c,'Area','Centroid','Orientation','Extent','Perimeter','ConvexArea','FilledArea','Solidity','Eccentricity','MajorAxisLength' ,'EquivDiameter','MinorAxisLength','MaxIntensity','MinIntensity','MeanIntensity','PixelList','PixelIdxList','SubArrayIdx','BoundingBox');
mindist = 2601;
cellStats = struct();
for i=1:numel(stats)
    if stats(i).Area > 3000 || stats(i).Area < 50
        continue;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ÜBERPRÜFEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif pdist([params.center;stats(i).Centroid]) < mindist %lastCentroid is doch falsch
        cellStats = stats(i);
        cellStats.bw_c = L;
        cellStats.I_c = I_c;
        mindist = pdist([params.center;stats(i).Centroid]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ÜBERPRÜFEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%     fprintf('hier soll keine 2601 oder so stehen:%i\n',cellStats.ConvexArea)
if numel(fields(cellStats)) == 0
    if usebw
        status = 3;
    else
        status = 4;
    end
elseif isnan(cellStats.Area)
    %     if cellStats.ConvexArea <= params.bmean+2*params.bstd
    if usebw
        cellStats.isbw = 1;
    else
        cellStats.isbw = 0;
    end
elseif cellStats.Area > 3000
    cellStats = struct();
    if usebw
        status = 1;
    else
        status = 3;
    end
elseif cellStats.Eccentricity > .99 && cellStats.Area > 3000
    status = 2;
    if usebw
        cellStats.isbw = 1;
    else
        cellStats.isbw = 0;
    end
else
    status = 1;
    if usebw
        cellStats.isbw = 1;
    else
        cellStats.isbw = 0;
    end
    %         cellStats = measureObjectAreaShape(I_c,L,cellStats);
end
end


