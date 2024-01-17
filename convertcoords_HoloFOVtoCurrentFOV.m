% xyintegration = updateSIrois_meso(hSI,holoRequest,clearprev, nstrips_orig, fullnpix_orig, fullxsize_orig, fullysize_orig, fullxcenter_orig, fullycenter_orig);
function xynew = convertcoords_HoloFOVtoCurrentFOV(hSI,xyorig, fullnpix_orig, fullxsize_orig, fullysize_orig, fullxcenter_orig, fullycenter_orig)

%%%%%%%%%%%% Convert holemasks to absolute um coords
[fullxpix_orig,fullypix_orig] = meshgrid(linspace(1,fullnpix_orig(1),fullnpix_orig(1)),...
    linspace(1,fullnpix_orig(2),fullnpix_orig(2)));
% fullxsize_orig = 500*nstrips_orig; %um
% fullysize_orig = 1000; %um
% fullxcenter_orig = 0; %um
% fullycenter_orig = -200; %um
[fullxum_orig,fullyum_orig] = meshgrid(...
    linspace(fullxcenter_orig-fullxsize_orig/2,fullxcenter_orig+fullxsize_orig/2,fullnpix_orig(1)),...
    linspace(fullycenter_orig-fullysize_orig/2,fullycenter_orig+fullysize_orig/2,fullnpix_orig(2))...
    );
fullxum_orig = round(fullxum_orig);
fullyum_orig = round(fullyum_orig);
xyum = NaN(size(xyorig));
for i=1:size(xyorig,1)
%     currind = fullxpix_orig==xyorig(i,1) & fullypix_orig==xyorig(i,2);
    currind = find(abs(fullxpix_orig-xyorig(i,1))==min(min(abs(fullxpix_orig-xyorig(i,1)))) ...
        & abs(fullypix_orig-xyorig(i,2))==min(min(abs(fullypix_orig-xyorig(i,2)))));
    xyum(i,:) = [fullxum_orig(currind),fullyum_orig(currind)];
end
%%%%%%%%%%%%% Holemasks converted to absolute um coordinates

%%%%%%%%%%%%% Convert holemasks to pixels in current/imaged sf/ROI
nstrips = length(hSI.hRoiManager.currentRoiGroup.rois);
nxpix = zeros(nstrips,1);
nypix = zeros(nstrips,1);
xsize = zeros(nstrips,1);
ysize = zeros(nstrips,1);
xcenter = zeros(nstrips,1);
ycenter = zeros(nstrips,1);
for n=1:nstrips
    currsf = hSI.hRoiManager.currentRoiGroup.rois(n).scanfields(1);
    nxpix(n) = currsf.pixelResolutionXY(1);
    nypix(n) = currsf.pixelResolutionXY(2);
    xsize(n) = currsf.sizeXY(1)*150; %um
    ysize(n) = currsf.sizeXY(2)*150; %um
    xcenter(n) = currsf.centerXY(1)*150; %um
    ycenter(n) = currsf.centerXY(2)*150; %um
end

fullnpix = [sum(nxpix),mean(nypix)];

[fullxpix,fullypix] = meshgrid(linspace(1,fullnpix(1),fullnpix(1)),...
    linspace(1,fullnpix(2),fullnpix(2)));
fullxsize = sum(xsize); %um
fullysize = mean(ysize); %um
fullxcenter = (xcenter(1)-xsize(1)/2) + ...
    abs((xcenter(1)-xsize(1)/2) - (xcenter(end)+xsize(end)/2))/2; %um
fullycenter = mean(ycenter); %um
[fullxum,fullyum] = meshgrid(...
    linspace(fullxcenter-fullxsize/2,fullxcenter+fullxsize/2,fullnpix(1)),...
    linspace(fullycenter-fullysize/2,fullycenter+fullysize/2,fullnpix(2))...
    );
fullxum = round(fullxum);
fullyum = round(fullyum);

xynew = NaN(size(xyorig));
for i=1:size(xyum,1)
%     currind = find(abs(fullxum-xyum(i,1))<=1 & abs(fullyum-xyum(i,2))<=1);
% changed by HS 220714 to make isequal(xyorig, xynew) true when ActualHoloFOV.roi is used
    currind = find(abs(fullxum-xyum(i,1))==min(min(abs(fullxum-xyum(i,1)))) & abs(fullyum-xyum(i,2))==min(min(abs(fullyum-xyum(i,2)))));
    [currindx,currindy] = ind2sub(size(fullxum),currind(1));
    %     currindx = abs(fullnpix(1)-currindx);
    xynew(i,:) = [fullxpix(currindx,currindy),fullypix(currindx,currindy)];
end
%%%%%%%%%%%%%% Holemasks converted to pixels in current ROI/FOV


