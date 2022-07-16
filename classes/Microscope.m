%% CGM_magicWandSegmentation package
% Class that defines the microscope used during the CGM experiments

% authors: Maelle Bénéfice, Guillaume Baffou
% affiliation: CNRS, Institut Fresnel
% date: Jul 31, 2019

classdef Microscope
    properties
        M
        dxSize
    end

    properties(Dependent)
        pxSize
    end

    methods

        function obj = Microscope(opt)
            arguments
                opt.M = 100
                opt.dxSize = 6.5e-6
            end
            obj.M=opt.M;
            obj.dxSize=opt.dxSize;

        end

        function val = get.pxSize(obj)
            val=obj.dxSize/obj.M;
        end
    end

end