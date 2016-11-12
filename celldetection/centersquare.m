function [bw_single_new,I_new,bw_single_new_resized,I_new_resized,Iorg_new,Iorg_new_resized] = centersquare(bw_single,bw,I,I_org,xs,xe,ys,ye,params)

maxbwsize = max(size(bw_single));

Imean = nanmean(I(~isinf(I)));
Iorgmean = nanmean(I_org(~isinf(I_org)));

zero1_bw = zeros(size(bw_single,1),maxbwsize);
zero1_I = (zero1_bw+1).*Imean;
zero1_Iorg = (zero1_bw+1).*Iorgmean;
bw_single = [zero1_bw bw_single zero1_bw];
bw = [zero1_bw bw zero1_bw];
I = [zero1_I I zero1_I];
I_org = [zero1_Iorg I_org zero1_Iorg];

zero2_bw = zeros(maxbwsize,size(bw_single,2));
zero2_I = (zero2_bw+1).*Imean;
zero2_Iorg = (zero2_bw+1).*Iorgmean;
bw_single = [zero2_bw; bw_single; zero2_bw];
bw = [zero2_bw; bw; zero2_bw];
I = [zero2_I; I; zero2_I];
I_org = [zero2_Iorg; I_org; zero2_Iorg];


assert(all(bw_single(:)==1| bw_single(:)==0),'binary image needed');
%first, figure out the centriod of the object

g00 = raw_moments(bw_single,0,0);
centroid = round([raw_moments(bw_single,1,0) raw_moments(bw_single,0,1) ]./g00 );

%find a circle around the oject
[X, Y] = meshgrid(1:size(bw_single,2),1:size(bw_single,1));
d = sqrt((X-centroid(1)).^2 +(Y-centroid(2)).^2) .* double(bw_single==1);
d = ceil(max(d(:))); %distance to the mfarthest point away from the centroid

%shrink the image to the box with size 2*d+1, this is where the circle fits
newX = (centroid(1)-d):(centroid(1)+d);
newY = (centroid(2)-d):(centroid(2)+d);
bw_single_new = logical(bw_single(newY,newX));
bw_new = logical(bw(newY,newX));
I_new = I(newY,newX);
Iorg_new = I_org(newY,newX);

newX_resized = (centroid(1)-floor(params.patchsize/2)):(centroid(1)+floor(params.patchsize/2));
newY_resized = (centroid(2)-floor(params.patchsize/2)):(centroid(2)+floor(params.patchsize/2));
bw_single_new_resized = logical(bw_single(newY_resized,newX_resized));
bw_new_resized = logical(bw(newY_resized,newX_resized));
I_new_resized = I(newY_resized,newX_resized);
Iorg_new_resized = I_org(newY_resized,newX_resized);

if params.normalize
    back = I_new(~bw_new);
    cell = I_new(bw_new);
    I_new = (I_new - mean(back(back~=0)))./std(cell(cell~=0));
    I_new = (255/2)+I_new.*(255/2);
    
    rback = I_new_resized(~bw_new_resized);
    rcell = I_new_resized(bw_new_resized);
    I_new_resized = (I_new_resized - mean(rback(rback~= 0)))./std(rcell(rcell~= 0));
    I_new_resized = uint8((255/2)+I_new_resized.*(255/2));

end



assert(raw_moments(bw_single_new,0,0)==g00,'something wrong g00 changed')
assert(size(bw_single_new,1)==size(bw_single_new,2),'image is not quadratic')


end

function outmom = raw_moments(im,i,j)
    outmom = sum(sum( ((1:size(im,1))'.^j * (1:size(im,2)).^i) .* im ));
end