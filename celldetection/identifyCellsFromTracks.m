function [Is_c_centered,Iorgs_c_centered,bws_c_centered, cellsizes, cellspeeds, label] = identifyCellsFromTracks(track,cellnr,params)

if params.doplot
    figure(1);
end

TPempty = [];
trackAttributes = fields(track);

cellidx = 1;

%initialize some base parameters
params.center = [fix(params.windowsize/2) fix(params.windowsize/2)];
params.lastCentroid = params.center;
params.currentCell = -1;
params.prevCell = -1;
params.prevSize = -1;

celltrack = struct();
for f =1:numel(trackAttributes)
    celltrack.(trackAttributes{f}) = track.(trackAttributes{f})(track.cellNr == cellnr);
end

%this is necessary as sometimes track timepoints are not sorted
tprange = min(celltrack.timepoint):max(celltrack.timepoint);



Is_c_centered = repmat({NaN},size(tprange));
Iorgs_c_centered = repmat({NaN},size(tprange));
bws_c_centered = repmat({NaN},size(tprange));
cellsizes = NaN(size(tprange));
absoluteX = NaN(size(tprange));
absoluteY = NaN(size(tprange));
idx = 1;

%iterate over all tracked timepoints
for t = tprange%(1:10:end)
    j = find(t == celltrack.timepoint);
    fprintf('\nProcessing cell %i - timepoint %i of %i\n',cellnr,find(t==tprange),size(tprange,2));
    
    params.windowsize = params.constwindowsize;
    params.center = [fix(params.windowsize/2) fix(params.windowsize/2)];
    
    tic
    
    oldx = celltrack.X(j);
    x = oldx;
    oldy = celltrack.Y(j);
    y = oldy;
    
    %load the image that shows the cell at the tracked timepoint and the
    %accompanying background
    I_org = im2double(imread(sprintf('%s%s/position%i/%s_pos%.4i_t%.5i_wl00.png',...
        params.datapath,params.experimentname,celltrack.positionIndex(j),params.experimentname,celltrack.positionIndex(j),celltrack.timepoint(j))));
    bg =imread(sprintf('%s%s/backgrounds/position%i/%s_pos%.4i_wl00_projbg.png',...
        params.datapath,params.experimentname,celltrack.positionIndex(j),params.experimentname,celltrack.positionIndex(j)));
    bg = double(bg)/(2^16-1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %Single Image Segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Normalize the image
    I = I_org./bg;
    I = I./mean(I(~isinf(I)));
    I_tp = I;
    
    %Empirically determined correction threshold that depends on the used
    %time lapse microscope.
    I(I>1.2) = 1.2;
    I(I<.7) = .7;
    
    %Try to segment the image. This might fail for various reasons so we
    %catch the error
    if params.segmentTotal
        try
            bw = segmentImage(I,params);
        catch me
            warning(me.message)
            idx = idx+1;
            continue
        end
        
        if params.dowatershed
            %Apply watershedding
            [~,bw,~] = watershedCell(I,bw);
            bw = logical(imfill(bw,'holes'));
        end
    else
        bw = imread(sprintf('%s%s/segmentation/position%i/%s_pos%.4i_t%.5i_wl00.png',...
            params.datapath,params.experimentname,celltrack.positionIndex(j),params.experimentname,celltrack.positionIndex(j),celltrack.timepoint(j)));
    end
    
    stats = regionprops(bw,I,'Eccentricity','Area','PixelIdxList','PixelValues');
    
    eccfilter = 0.99;
    ecccombfilter = 0.7;
    maxsizecombfilter = 700;
    minsize = 40;
    maxsize = 1000;
    
    tofilt=([stats.Eccentricity]>eccfilter | [stats.Area]<minsize | [stats.Area]>maxsize | ([stats.Eccentricity]>ecccombfilter & [stats.Area]>maxsizecombfilter)); %| ([stats.Eccentricity]<0.5 & [stats.Area]>maxSize);
    for id = find(tofilt)
        bw([stats(id).PixelIdxList]) =0;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %Single Cell Segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    [I_c,xborder,yborder,xcut,ycut] = cropImage(I_tp,x,y,params.windowsize);
    params.relativetrackPoint = [params.center(1)+xcut,params.center(2)+ycut];
    [xs,ys] = size(I_c);
    if isempty(I_c) || (xs*ys) <= 50
        warning('Empty cropped image, skipping timepoint!')
        idx = idx+1;
        continue
    end
    
    bw = logical(imfill(bw,'holes'));
    bw_c = cropImage(bw,x,y,params.windowsize);
    
    while (sum(sum(bw_c))) == 0 && params.windowsize <= params.maxwindowsize
        params.windowsize = params.windowsize + params.growwindow;
        params.center = [fix(params.windowsize/2) fix(params.windowsize/2)];
        
        [I_c,xborder,yborder] = cropImage(I,x,y,params.windowsize);
        bw_c = cropImage(bw,x,y,params.windowsize);
    end
    
    if params.windowsize > params.maxwindowsize
        fprintf('Track point is very far off or cell was not segmented., skipping this timepoint!\n');
        idx = idx+1;
        continue;
    elseif params.windowsize > 4*params.constwindowsize
        fprintf('Track point is very far off or cell was not segmented., skipping this timepoint!\n');
        if doplot
            clf
            subplot(2,3,[1 2 4 5])
            imagesc(I)
            colormap gray
            hold on
            plot(oldx,oldy,'b+')
            plot(celltrack.X(j),celltrack.Y(j),'r+')
            drawnow;
        end
        idx = idx+1;
        continue;
    end
    
    while sum(sum(imclearborder(bw_c))) == 0 && params.windowsize <= params.maxwindowsize
        params.windowsize = params.windowsize + params.growwindow;
        %blob is touching the border, extend the image size
        tempprops = regionprops(bw_c,'Centroid');
        mindist = inf;
        for object=1:numel(tempprops)
            if pdist([params.center;tempprops(object).Centroid]) < mindist
                tempprop = tempprops(object);
                
                mindist = pdist([params.lastCentroid;tempprops(object).Centroid]);
            end
        end
        
        tempoffsetX = round(tempprop.Centroid(1)-params.windowsize/2);
        tempoffsetY = round(tempprop.Centroid(2)-params.windowsize/2);
        imagesize = size(I); % X und Y sind rumgedreht!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        %set the X and Y for one position
        if x+tempoffsetX > imagesize(2)
            x = imagesize(1);
        elseif x+tempoffsetX < 0
            x = 0;
        else
            x = x+tempoffsetX;
        end
        if y+tempoffsetY > imagesize(1)
            y = imagesize(2);
        elseif y+tempoffsetY < 0
            y = 0;
        else
            y = y+tempoffsetY;
        end
        %set the X and Y for the whole plate should not happen that it is bigger or smaller than possible?!?
        I_c = cropImage(I,x,y,params.windowsize);
        [bw_c,xborder,yborder,~,~] = cropImage(bw,x,y,params.windowsize);
        celltrack.X(j) = x;
        celltrack.Y(j) = y;
        
    end
    
    if sum(sum(bw_c)) == 0
        disp('No segmented cell found in bw patch')
        idx = idx+1;
        continue;
    else
        if params.normalize
            back = I_c(~bw_c);
            cell = I_c(bw_c);
            I_c = (I_c - mean(back(:)))./std(cell(:));
            %                 I_c = (I_c - mean(back(:)));
            I_c = uint8((255/2)+I_c.*(255/2));
            I_c(I_c < 0) = 0;
            I_c(I_c > 255) = 255;
        end
        [props,params,~] = segmentSingleTimepoint(I_c,bw_c,params);
        
        
        % only continue if a blob was identified
        if numel(fields(props)) == 0
            disp('No segmented cell found in bw patch')
            idx = idx+1;
            continue;
        else
            
            %if find(t==tprange) == 51
            %    'h'
            %end
            
            %reconstruct the found object on the big image
            xoffset = round((params.windowsize-size(bw_c,2))/2)*xborder;
            yoffset = round((params.windowsize-size(bw_c,1))/2)*yborder;
            props.CroppedCentroid = props.Centroid;
            props.Centroid(1) = round(celltrack.X(j)-(size(bw_c,2)/2-props.CroppedCentroid(1))+xoffset)-1;
            props.Centroid(2) = round(celltrack.Y(j)-(size(bw_c,1)/2-props.CroppedCentroid(2))+yoffset)-1;
            
            props.CroppedBoundingBox = props.BoundingBox;
            props.BoundingBox(1) = fix(celltrack.X(j)-(size(bw_c,2)/2-props.CroppedBoundingBox(1))+xoffset);
            props.BoundingBox(2) = fix(celltrack.Y(j)-(size(bw_c,1)/2-props.CroppedBoundingBox(2))+yoffset);
            
            props.CroppedPixelList = props.PixelList;
            props.PixelList(:,1) = round(props.CroppedPixelList(:,1)+xoffset+(celltrack.X(j)-size(bw_c,2)/2))-1;
            props.PixelList(:,2) = round(props.CroppedPixelList(:,2)+yoffset+(celltrack.Y(j)-size(bw_c,1)/2))-1;
            
            %build new big bw image
            bw_single = false(size(bw));
            bw_single(sub2ind(size(bw),props.PixelList(:,2),props.PixelList(:,1))) = true;
            
            %build centered and cleaned cropped images
            blob_bb = props.BoundingBox;
            
            props.bw = bw_single;
            xs = max(blob_bb(2),1);
            xe = blob_bb(2)+min(blob_bb(3),blob_bb(4));
            ys = max(blob_bb(1),1);
            ye = blob_bb(1)+min(blob_bb(3),blob_bb(4));
                        
            [bw_c_centered, I_c_centered, bw_c_centered_resized, I_c_centered_resized, Iorg_c_centered,Iorg_c_centered_resized] = centersquare(bw_single,bw,I,I_org,xs,xe,ys,ye,params);
        end
        
        Is_c_centered{idx} = I_c_centered_resized;
        Iorgs_c_centered{idx} = Iorg_c_centered_resized;
        bws_c_centered{idx} = bw_c_centered_resized;
        cellsizes(idx) = props.Area;
        
        [absXcorr,absYcorr] = computeAbsoluteCoordinates(oldx, oldy, props.Centroid, celltrack.absoluteX(j), celltrack.absoluteY(j), params.mperp);
        absoluteX(idx) = absXcorr;
        absoluteY(idx) = absYcorr;
        
        % plot the stuff if wanted
        if params.doplot
            clf
            subplot(2,3,[1 2 4 5])
            imagesc(I)
            colormap gray
            hold on
            plot(props.Centroid(1),props.Centroid(2),'g+')
            plot(oldx,oldy,'b+')
            plot(celltrack.X(j),celltrack.Y(j),'r+')
            contour(bw_single,[1 1],'Color','g')
            bb = props.BoundingBox;
            bb(1) = bb(1)-5;
            bb(2) = bb(2)-5;
            bb(3) = bb(3)+5;
            bb(4) = bb(4)+5;
            rectangle('Position',bb,'EdgeColor','b')
            rectangle('Position',[celltrack.X(j)-(params.windowsize/2) ...
                celltrack.Y(j)-(params.windowsize/2) params.windowsize ...
                params.windowsize],'EdgeColor','w')
            
            subplot(2,3,3)
            imagesc(I_c_centered_resized)
            caxis([min(I_c_centered_resized(:)) max(I_c_centered_resized(:))])
            hold on
            contour(bw_c_centered_resized,[1 1],'g')
            
            axis equal
            hold off
            subplot(2,3,6)
            
            plot(cellsizes,'.','LineWidth',5)
            xlim([1 numel(tprange)])
            
            ylim([0 400]);
            xlabel('Timepoint')
            ylabel(['Cell Size'])
            title(celltrack.positionIndex(j));
            drawnow
        end
    end
    clear props
    idx = idx+1;
end
cellspeeds = computeCellSpeed(absoluteX,absoluteY, celltrack.absoluteTime);

type = unique(celltrack.type);
invtype = unique(celltrack.invtype_alt);
if type == 1
        label = 'MEP';
elseif type == 22
        label = 'GMP';
elseif invtype == 1
    label = 'MEP';
elseif invtype == 2
    label = 'GMP';
end
end

function [absXcorr,absYcorr] = computeAbsoluteCoordinates(oldx, oldy, centroid, absoluteX, absoluteY, mperp)
    
    offsetX = oldx - centroid(1);
    offsetY = oldy - centroid(2);
    
    absXcorr = absoluteX + (offsetX * mperp);
    absYcorr = absoluteY +(offsetY * mperp);

end

function movementvect = computeCellSpeed(X,Y,abstime)
movementvect = nan(numel(X),1);
movementvect(1) = 0;
for i=1:numel(X)-1
    movepoint =  sqrt(((X(i)-X(i+1))^2)+(Y(i)-Y(i+1))^2);
    movepoint = movepoint/(abstime(i+1)-abstime(i));
    movementvect(i+1) = movepoint;
end
end
