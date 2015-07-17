classdef ECPUHF < SCF.UHF & SCF.ECPRHF
    
    methods
        
        function obj = ECPUHF(properties)
            obj = obj@SCF.ECPRHF(properties);
            obj = obj@SCF.UHF(properties);
        end
        
        function [densVec, orbital] = HarrisGuess(obj)
            orbital = {SCF.ECPRHF.G09ReadMatrix('HarrisGuessMOAlpha'), ...
                SCF.ECPRHF.G09ReadMatrix('HarrisGuessMOBeta')};
            order = SCF.ECPRHF.G09ToPsi4BasisOrder(obj.matpsi2.BasisSet_ShellNumFunctions());
            orbital{1} = orbital{1}(order, :);
            orbital{2} = orbital{2}(order, :);
            densVec = obj.OrbToDensVec(orbital);
        end
        
    end
    
end