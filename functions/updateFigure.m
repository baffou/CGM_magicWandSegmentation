%% CGM_magicWandSegmentation package
% In case a list of images is loaded, this function updates the T and OPD
% images of the IM.figure() according to the property hfig.UserData.NumIm.

% authors: Maelle Bénéfice, Guillaume Baffou
% affiliation: CNRS, Institut Fresnel
% date: Jan 18, 2021

function updateFigure(hfig)
    if strcmp(hfig.UserData.images,'both')
        hfig.UserData.ax(1).Children.CData=hfig.UserData.IM(hfig.UserData.NumIm).T;
        hfig.UserData.ax(2).Children.CData=hfig.UserData.IM(hfig.UserData.NumIm).OPD;
        drawnow
        title(hfig.UserData.ax(1),['transmission ' num2str(hfig.UserData.NumIm) '/' num2str(numel(hfig.UserData.IM))])
        title(hfig.UserData.ax(2),['wavefront ' num2str(hfig.UserData.NumIm) '/' num2str(numel(hfig.UserData.IM))])
    elseif strcmp(hfig.UserData.images,'T')
        hfig.UserData.ax(1).Children.CData=hfig.UserData.IM(hfig.UserData.NumIm).T;
        drawnow
        title(hfig.UserData.ax(1),['transmission ' num2str(hfig.UserData.NumIm) '/' num2str(numel(hfig.UserData.IM))])
    elseif strcmp(hfig.UserData.images,'OPD')
        hfig.UserData.ax(1).Children.CData=hfig.UserData.IM(hfig.UserData.NumIm).OPD;
        drawnow
        title(hfig.UserData.ax(1),['wavefront ' num2str(hfig.UserData.NumIm) '/' num2str(numel(hfig.UserData.IM))])
    end
    drawnow
end
