classdef ECPUHF < SCF.UHF & SCF.ECPRHF
    
    methods
        
        function obj = ECPUHF(properties)
            obj = obj@SCF.ECPRHF(properties);
            obj = obj@SCF.UHF(properties);
        end
        
        function [densVec, orbital] = HarrisGuess(obj)
            orbital = {SCF.ECPRHF.G09ReadMatrix('HarrisGuessMOAlpha'), ...
                SCF.ECPRHF.G09ReadMatrix('HarrisGuessMOBeta')};
            densVec = obj.OrbToDensVec(orbital);
        end
        
    end
    
end