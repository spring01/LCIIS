classdef ECPRKS < SCF.RKS & SCF.ECPRHF
    
    methods
        
        function obj = ECPRKS(properties, dft)
            obj = obj@SCF.ECPRHF(properties);
            obj = obj@SCF.RKS(properties, dft);
        end
        
    end
    
end