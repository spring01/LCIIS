classdef G09USCF < SCF.G09RSCF & SCF.UHF
    
    methods
        
        function obj = G09USCF(info)
            obj = obj@SCF.G09RSCF(info);
            matrices = SCF.G09RSCF.G09ReadMatrices({'overlap', 'coreHamilt'});
            properties.overlapMat = matrices{1};
            properties.coreHamilt = matrices{2};
            scalars = SCF.G09RSCF.G09ReadScalars({'NumElectrons', 'NucRepEnergy'});
            properties.numElectrons = scalars{1};
            properties.nucRepEnergy = scalars{2};
            properties.matpsi2 = [];
            obj@SCF.UHF(properties);
            obj.info = info;
        end
        
        function [densVec, orbital] = HarrisGuess(obj)
            info_ = obj.info;
            info_.harris = [];
            SCF.G09RSCF.RunG09(info_);
            matrices = SCF.G09RSCF.G09ReadMatrices({'HarrisGuessMOAlpha', 'HarrisGuessMOBeta'});
            orbital = matrices;
            densVec = obj.OrbToDensVec(orbital);
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(obj, orbital)
            obj.info.orbAlpha = orbital{1};
            obj.info.orbBeta = orbital{2};
            SCF.G09RSCF.RunG09(obj.info);
            matrices = SCF.G09RSCF.G09ReadMatrices({'fockAlpha', 'fockBeta'});
            fockVec = [reshape(matrices{1}, [], 1), reshape(matrices{2}, [], 1)];
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            energy = obj.SCFEnergy@SCF.G09RSCF(fockVec, densVec);
        end
        
    end
    
end
