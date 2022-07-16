%% CGM_magicWandSegmentation package
% function that deals with the magic wand on a CGM image

% authors: Maelle Bénéfice, Guillaume Baffou
% affiliation: CNRS, Institut Fresnel
% date: Jan 18, 2021

function OV = magicWand_OV(hfig)

Nim=numel(hfig.UserData.IM);
OV=zeros(Nim,1);
for ii=1:Nim
    hfig.UserData.NumIm=ii;
    updateFigure(hfig)

    IM = hfig.UserData.IM(hfig.UserData.NumIm);
    
    bkTh = 3;
    step = 5;
    nmax=20;
    
    OPD = IM.OPD;

    YLim = hfig.UserData.ax(1).XLim;
    XLim = hfig.UserData.ax(1).YLim;
    
    fail = 1;
    while fail==1 % until space bar or c are pressed
        OPDcrop = OPD;
        [mask,maskRemove,fail,xList,yList] = magicWand_scrollbar(hfig);
    end %end while fail==1
    
    maskRef = mask;
    
    %% OV calculation
    
    %when mask smaller than bacteria
    taillePx = IM(1).Microscope.pxSize;
    N = round(nmax/2);
    alpha0 = zeros(1,N+nmax);
    OV0 = zeros(1,N+nmax);
    xcoord = zeros(1,N+nmax);
    maskList = cell(N+nmax,1);
    ringBList = cell(N+nmax,1);
    for i = 1:N  % loop over mask size from smallest to size of ref.
        Wb = wiener2(mask,[2*bkTh+1,2*bkTh+1]); %sets the width of the background
        Ws = wiener2(mask,[2*step+1,2*step+1]);       %sets the step for the calculations
        maskb = double((Wb>0).*maskRemove);
        masks = double((Ws==1).*maskRemove);
        ringB = maskb-mask;
        
        backring = ringB.*OPDcrop;
        background = sum(backring(:))/(sum(ringB(:))); %mean backgound over the ring
        OPDn = OPDcrop-background;
        
        OV0(N-i+1) = sum(OPDn(:).*mask(:))*taillePx.*taillePx;
        xcoord(N-i+1) = (N-i+1)/N;
        maskList{N-i+1} = mask;
        ringBList{N-i+1} = ringB;
        mask = masks;
    end
    
    mask = maskRef;
    for i = N:N+nmax  % loop over mask size when mask bigger than the ref mask
        Wb = wiener2(mask,[2*bkTh+1,2*bkTh+1]);
        Ws = wiener2(mask,[2*step+1,2*step+1]);
        maskb = double((Wb>0).*maskRemove);
        masks = double((Ws>0).*maskRemove);
        ringB = maskb-mask;
        
        backring = ringB.*OPDcrop;
        background = sum(backring(:))/sum(ringB(:));
        OPDn = OPDcrop - background;
        
        
        OV0(i) = sum(OPDn(:).*mask(:))*taillePx.*taillePx;
        xcoord(i) = (i)/N;
        maskList{i} = mask;
        ringBList{i} = ringB;
        mask = masks;
    end
    
    hfig2 = figure;
    hfig2.UserData.goon = 0;
    hfig2.UserData.p1 = [];
    fullwidth
    ha1 = subplot(1,2,1);
    hold on
    plot(xcoord,real(alpha0))
    plot(xcoord,imag(alpha0))
    plot(xcoord,OV0)
    legend({'real(\alpha)','imag(\alpha)','OV'})
    ha2 = subplot(1,2,2);
    A = zeros([size(OPDn) 3]);
    
     maxval = max(max(imgaussfilt(OPDn,10)));
     minval = min(min(imgaussfilt(OPDn,10)));
    
    A(:,:,1) = (OPDn-minval)/(maxval-minval);
    A(:,:,2) = (OPDn-minval)/(maxval-minval);
    A(:,:,3) = (OPDn-minval)/(maxval-minval);
    
    A(:,:,1) = (mask~=0); % pixel corresponding to the mask set to 1 in red channel
    image(A)
    set(gca,'YDir','Reverse')
    drawnow
    hfig2.UserData.ha1 = ha1;
    hfig2.UserData.ha2 = ha2;
    set(ha1,'ButtonDownFcn',@(src,evt)setp1p2(ha1));
    set(ha1.Children,'PickableParts','none');
    %set(hfig2,'WindowButtonMotionFcn','disp(get (gca, ''CurrentPoint''))')
    set(hfig2,'WindowButtonMotionFcn',@(src,evt)LiveMaskDisplay(ha1,ha2,maskList,N))
    %[xp,~] = ginput(2);
    pause(0.1)
    
    
    while(hfig2.UserData.goon==0)
    pause(0.1)
    end
    
    xp(1) = hfig2.UserData.p1;
    xp(2) = hfig2.UserData.p2;
    close(hfig2)
    
    pxmin = round(min(N*xp));
    pxmax = min(round(max(N*xp)),length(alpha0)); % in case the user clicks too far in x
    if pxmax<pxmin
        error('The second click must correspond to a higher x value.')
    end
    maskMeas = maskList{round(mean([pxmin pxmax]))};
    ringBMeas = ringBList{round(mean([pxmin pxmax]))};
    
    OV(ii) = mean(OV0(pxmin:pxmax));
    
    fprintf('OV:\t%.4g pg\n',OV(ii)*1e18);
    
    folder = 'results';
    if ~isfolder(folder)
        mkdir(folder)
    end
    
    %% save the images
    if ~isfolder([folder '/masks'])
        mkdir([folder '/masks'])
    end
    %% Save data
    
    writematrix(OPDn,[folder '/masks/OPDcrop_' IM.OPDfileName])
    writematrix(maskMeas,[folder '/masks/Mask_' IM.OPDfileName])
    writematrix(xcoord,[folder '/masks/pcoord_' IM.OPDfileName])
    writematrix(OV0,[folder '/masks/OVprofile_' IM.OPDfileName])
    minOPD = min(min(imgaussfilt(OPDn,2)));
    maxOPD = max(max(imgaussfilt(OPDn,2)));
    imwrite((OPDn-minOPD)*255/(maxOPD-minOPD),phase1024(256),[folder '/masks/OPDcrop_' IM.OPDfileName(1:end-4) '.tif'])
    imwrite(maskMeas*255,gray(256),[folder '/masks/Mask_' IM.OPDfileName(1:end-4) '.tif'])
    imwrite(ringBMeas*255,gray(256),[folder '/masks/RingB_' IM.OPDfileName(1:end-4) '.tif'])
    
    figure(hfig)
    
    %% Save csv file

        %% save the csv file
    % if csv file does not already created, write the first line with the row's names.
    if ~exist([folder '/data.csv'],'file') 
        Tab0 = table({'file name','time','OV (pg)','XY'});
        writetable(Tab0,[folder '/data.csv'],'WriteVariableNames',0);
    end
    efactor=1e-18;
    % defines the line to write
    posList = [XLim(1)+xList,YLim(1)+yList].';
    date = generateDatedFileName();
    Tab = table({IM.OPDfileName},{date},{sprintf('%.4g',OV(ii)/efactor)},{num2str(posList(:).')});
    % writes this line in a temporary csv file
    %fileNamecsv = [folder '/data_' date '.csv'];
    writetable(Tab,[folder '/data.csv'],'WriteVariableNames',0,'WriteMode','Append');
    % concatenates the current data.csv file with the new one 
    %copyfile([folder '/data.csv'],[folder '/data_ref.csv'])
    %command = ['cat ' folder '/data_ref.csv ' fileNamecsv ' > ' folder '/data.csv'];
    %system(command);
    %delete(fileNamecsv)
    delete([folder '/data_ref.csv'])


end

end



