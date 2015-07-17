classdef ECPUKS < SCF.UKS & SCF.ECPUHF
    
    methods
        
        function obj = ECPUKS(properties, dft)
            obj = obj@SCF.ECPUHF(properties);
            obj = obj@SCF.UKS(properties, dft);
        end
        
    end
    
end