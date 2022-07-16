addpath(genpath(pwd))
clear
close all

%% import images
MI = Microscope('M',60,'dxSize',6.5e-6);

IM = ImageQLSI(4);
for ii = 1:4
    TfileName = ['data/T' textk(ii) '.txt'];
    OPDfileName = ['data/OPD' textk(ii) '.txt'];
    IM(ii) = ImageQLSI(TfileName,OPDfileName,MI);
end

IM.figure('images','OPD')
IM.figure('axes','um','images','T')
IM.figure()

IMc = IM.crop("Center","Manual");


imwrite(IM(1).OPD,"OPD.tif")

%% process images


close all
h = IMc.figure('axes','um');
OV = magicWand_OV(h);



