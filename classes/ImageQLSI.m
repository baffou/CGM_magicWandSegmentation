%% CGM_magicWandSegmentation package
% Class that defines an experimental image acquired by CGM

% authors: Maelle Bénéfice, Guillaume Baffou
% affiliation: CNRS, Institut Fresnel
% date: Jul 31, 2019

classdef ImageQLSI
    
%    properties(GetAccess = public, SetAccess = private)
    properties(Access = public, Hidden)
        T       % matrix of path
        OPD     
        Microscope Microscope
    end
    
    properties(GetAccess = public, SetAccess = private)
        TfileName
        OPDfileName
    end
    
    properties(Dependent)
        Nx
        Ny
    end
    
    properties(Hidden,SetAccess = private)
        Fcrops %3-cell FcropParameters object crops{1:3}
        processingSoftware
    end
    
    methods
        function obj = ImageQLSI(INT,OPD,MI)
                
            if nargin==1 %ImageQLSI(n)
                if isnumeric(INT)
                    n=INT;
                    obj = repmat(ImageQLSI(),n,1);
                else
                    error('wrong input type')
                end
            elseif nargin==0 %ImageQLSI()
            elseif nargin==3 %ImageQLSI(T,OPD,MI)
                if isnumeric(INT) && isnumeric(OPD) % if both are matrices, they should be of the same size
                    if (sum(size(INT)==size(OPD))~=2)
                        size(INT)
                        size(OPD)
                        error('OPD and T images must have the same size.')
                    end
                end

                if ischar(INT)
                    obj.T = readmatrix(INT);
                    obj.TfileName=extractFileName(INT);
                elseif isnumeric(INT)
                    obj.T=INT;
                else
                    error('wrong input')
                end

                if ischar(OPD)
                    obj.OPD = readmatrix(OPD);
                    obj.OPDfileName=extractFileName(OPD);
                elseif isnumeric(OPD)
                    obj.OPD=OPD;
                else
                    error('wrong input')
                end

                obj.Microscope = MI;
                
            else
                error('Not the proper number of inputs')
            end
        end

        function val = get.Nx(obj)
            val = size(obj.T,2);
        end

        function val = get.Ny(obj)
            val = size(obj.T,1);
        end
        
        function obj = OPDhighPass(obj,nCrop)
            if nargin==1
                nCrop = 3;
            end
            %imageT_F=fftshift(fft2(obj.T));
            imageOPD_F = fftshift(fft2(obj.OPD));
            
            ix0 = floor(obj.Nx/2)+1;
            iy0 = floor(obj.Ny/2)+1;
            %imageT_F(iy0-nCrop:iy0+nCrop,ix0-nCrop:ix0+nCrop) = 0;
            imageOPD_F(iy0-nCrop:iy0+nCrop,ix0-nCrop:ix0+nCrop) = 0;
            %imageT2 = real(ifft2(ifftshift(imageT_F)));
            imageOPD2 = real(ifft2(ifftshift(imageOPD_F)));
            %obj.T = 1+imageT2;
            obj.OPD = imageOPD2;
            
            
        end
        
        function obj = ZernikeRemove(obj,n,m,r0)
            if nargin==1
                n = 1;
                m = 1;
            end
            if nargin==2 % if only n specified, remove all the Zernike polynoms up to n
                
                for ii = 1:n
                    for jj = 1:ii
                        m = mod(ii,2)+2*(jj-1);
                        obj.OPD = ZernikeRemoval(obj.OPD,ii,m);
                    end
                end
                
            end
            if nargin==4
                obj.OPD = ZernikeRemoval(obj.OPD,n,m,r0);
            else
                obj.OPD = ZernikeRemoval(obj.OPD,n,m);
            end
            %obj.OPD = ZernikeRemoval(obj.OPD,n,m);
            
        end
        
        function val = DM(obj)
            % returns the dry mass density in pg/µm^2
            val=5.56*obj.OPD;
        end

        function obj = crop(obj,opt)
            arguments
                obj
                opt.xy1 = 0
                opt.xy2 = 0
                opt.Center {mustBeMember(opt.Center,{'Auto','Manual','auto','manual'})} = 'Auto'
                opt.Size = 'Manual'
            end
            if opt.xy1~=0 && opt.xy2~=0  % crop('xy1', [a,b], 'xy2', [a,b])
                if numel(opt.xy1)~=2
                    error('First input must be a 2-vector (x,y)')
                end
                if numel(opt.xy2)~=2
                    error('Second input must be a 2-vector (x,y)')
                end
                x1 = min([opt.xy1(1),opt.xy2(1)]);
                x2 = max([opt.xy1(1),opt.xy2(1)]);
                y1 = min([opt.xy1(2),opt.xy2(2)]);
                y2 = max([opt.xy1(2),opt.xy2(2)]);
            else
                if strcmpi(opt.Center,'Manual') % crop('Center','Manual')
                    h = obj(1).figure(images='OPD');
                    [xc, yc] = ginput(1);
                    close(h)
                else
                    xc = obj(1).Nx/2;
                    yc = obj(1).Ny/2;
                end
                if strcmpi(opt.Size,'Manual') % crop('Center','Manual')
                    h = obj(1).figure(images='OPD');
                    rr = rectangle();
                    set (h, 'WindowButtonMotionFcn', @(src,event)drawRect(h,rr,xc,yc))
                    set (h, 'WindowButtonDownFcn', @(src,event)getPoint(h))
                    h.UserData.getPoint=0;
                    while h.UserData.getPoint==0
                    pause(0.1)
                    end
                    Size = h.UserData.side;
                    close (h)
                    x1 = floor(xc-Size/2+1);
                    x2 = floor(xc-Size/2+Size);
                    y1 = floor(yc-Size/2+1);
                    y2 = floor(yc-Size/2+Size);
                elseif isnumeric(opt.Size)
                    x1 = floor(xc-opt.Size/2+1);
                    x2 = floor(xc-opt.Size/2+opt.Size);
                    y1 = floor(yc-opt.Size/2+1);
                    y2 = floor(yc-opt.Size/2+opt.Size);
                end

            end
            
            for io = 1:numel(obj)
                temp=obj(io).T(y1:y2,x1:x2); % temp variable to avoid importing the matrix twice for the calculation of Nx and Ny when it is stored in a file.
                obj(io).T   = temp; 
                obj(io).OPD = obj(io).OPD(y1:y2,x1:x2);
            end

            function drawRect(h,rr,xc,yc)
                C = get (gca, 'CurrentPoint');
                side = floor(max(abs([2*(C(1)-xc), 2*(C(3)-yc)])));
                set(rr,'Position',[xc-side/2 yc-side/2 side side])
                h.UserData.side = side;
            end
            
            function getPoint(h)
                h.UserData.getPoint = 1;
            end
            
        end
          
        function h=figure(obj,opt)
            arguments
                obj
                opt.axes = 'px'
                opt.images  char {mustBeMember(opt.images,{'both','T','OPD'})} = 'both'
            end
            
            if strcmpi(opt.axes,'px')
                xx=1:obj(1).Nx;
                yy=1:obj(1).Ny;
            elseif strcmpi(opt.axes,'um')
                xx=(0:obj(1).Nx-1)*obj(1).Microscope.pxSize*1e6;
                yy=(0:obj(1).Ny-1)*obj(1).Microscope.pxSize*1e6;
            else
                error('Not a proper axes input')
            end

            % figure plot
            hfig = figure;
            hfig.UserData.NumIm = 1;
            hfig.UserData.images=opt.images;
            hfig.WindowButtonDownFcn=@clickFigureEvent;
            hfig.UserData.IM = obj;

            if strcmpi(opt.images,'T') || strcmpi(opt.images,'OPD')
                Nim=1;
            else
                Nim=2;
            end
%            hfig.UserData.ax=zeros(Nim,1);
            im=0;
            if strcmpi(opt.images,'T') || strcmpi(opt.images,'both')
                im=im+1;
                hfig.UserData.ax(im)=subplot(1,Nim,im);
                imagesc(xx,yy,obj(1).T)
                colormap(gca,'Gray')
                set(gca,'DataAspectRatio',[1,1,1])
                title('transmission')
                colorbar('southoutside');
                xlabel(opt.axes)
                ylabel(opt.axes)
                set(gca,'YDir','normal')

            end

            if strcmpi(opt.images,'OPD') || strcmpi(opt.images,'both')
                im=im+1;
                hfig.UserData.ax(im)=subplot(1,Nim,im);
                imagesc(xx,yy,obj(1).OPD*1e9)
                colormap(gca,phase1024)
                set(gca,'DataAspectRatio',[1,1,1])
                set(gca, 'YDir','reverse')
                title(['wavefront 1/' num2str(numel(hfig.UserData.IM))])
                c=colorbar('southoutside');
                c.Label.String = '(nm)';
                xlabel(opt.axes)
                ylabel(opt.axes)
                set(gca,'YDir','normal')    
                posFig=get(gcf,'Position');
            end

            if Nim==2
                set(gcf,'Position',[posFig(1)/2 posFig(2) 2*posFig(3) posFig(4)])
            end

            hfig.UserData.roiIN{1} = [];
            hfig.UserData.roiOUT{1} = [];

            if nargout
                h=gcf;
            end

            function clickFigureEvent(hfig)
                button=get(gcf,'SelectionType');
                if strcmp(button,'normal')
                    hfig.UserData.NumIm=hfig.UserData.NumIm+1;
                elseif strcmp(button,'alt')
                    hfig.UserData.NumIm=hfig.UserData.NumIm-1;
                end
                if hfig.UserData.NumIm<1
                    hfig.UserData.NumIm=1;
                    warning('You reach image #1, you cannot go below')
                elseif  hfig.UserData.NumIm>numel(hfig.UserData.IM)
                    hfig.UserData.NumIm=numel(hfig.UserData.IM);
                    warning(['There is no more than ' num2str(numel(hfig.UserData.IM)) ' images'])

                else
                    updateFigure(hfig)
                end

                    
            end


        end

    end
    
    
end








