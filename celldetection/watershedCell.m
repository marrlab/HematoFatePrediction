function [L,bw2,bw_woMerge] = watershedCell(I,bw)

warning('OFF','images:imfindcircles:warnForLargeRadiusRange')
warning('OFF','images:imfindcircles:warnForSmallRadius')

Rmin = 3;
Rmax = 25;

[centers,~] = imfindcircles(I,[Rmin Rmax],'ObjectPolarity','dark');
MaximaImage = false(size(bw));
for idx=1:size(centers,1)
    MaximaImage(round(centers(idx,2)),round(centers(idx,1))) = true;
end
MaximaImage = imdilate(MaximaImage,strel('disk',1));

if ~exist('DistanceTransformedImage','var')
    DistanceTransformedImage = bwdist(~bw);
end
Overlaid = imimposemin(-DistanceTransformedImage,MaximaImage);

wb = watershed(Overlaid) > 0;
bw = bw.*wb;

bw = bwlabel(bw);

bw_woMerge = bwlabel(bw > 0);
bw_woMerge(bw_woMerge > 0) = 1;


bw2 = bw_woMerge;
L = bw2;

warning('ON','images:imfindcircles:warnForLargeRadiusRange')
warning('ON','images:imfindcircles:warnForSmallRadius')
end
