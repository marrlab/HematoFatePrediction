function bw = segmentImage(I,params)

I = I-min(min(I));
I = I./max(max(I));
I = im2uint8(I);
                
msers = linearMser(im2uint8(I),params.delta,params.minsize,params.maxsize,1,true);
bw = zeros(size(I));
for i=1:numel(msers)
    bw(msers{i}) = 1;
end

bw = imclose(bw,strel('disk',2));
bw = bwmorph(bw,'bridge');
bw = logical(imfill(bw,'holes'));

bw = imopen(bw,strel('disk',2));

end
