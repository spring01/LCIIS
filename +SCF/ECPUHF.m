classdef ECPUHF < SCF.UHF & SCF.ECPRHF
    
    methods
        
        function obj = ECPUHF(properties)
            obj = obj@SCF.ECPRHF(properties);
            obj = obj@SCF.UHF(properties);
        end
        
    end
    
end