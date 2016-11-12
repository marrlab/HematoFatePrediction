function [I_c,xborder,yborder,xcut,ycut] = cropImage(image,x,y,windowsize)
    pluspixels = windowsize/2;
    minuspixels = windowsize/2;
    xborder = 1;
    yborder = 1;
    xcut = 0;
    ycut = 0;
    %try to determine if x or y are a t border
    if x+pluspixels <= size(image,2)
        xp = x+pluspixels-1;
    else
        xp = size(image,2);
        xborder = -1; 
        xcut = (x+pluspixels-1)-size(image,2);
    end
    if x-minuspixels >= 1
        xn = x-minuspixels;
    else
        xn = 1;
        xborder = 1;
        xcut = x-minuspixels;
    end
    if y+pluspixels <= size(image,1)
        yp = y+pluspixels-1;
    else
        yp = size(image,1);
        yborder = -1;
        ycut = (y+pluspixels-1)-size(image,1);
    end
    if y-minuspixels >= 1
        yn = y-minuspixels;
    else
        yn = 1;
        yborder = 1;
        ycut = y-minuspixels;
    end

    I_c = image(yn:yp,xn:xp); %matlab has switched x and y labels... 
end