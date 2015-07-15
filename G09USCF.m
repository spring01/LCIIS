classdef G09USCF < G09RSCF & UHF
    
    methods
        
        function obj = G09USCF(info)
            obj = obj@G09RSCF(info);
            matrices = G09RSCF.G09ReadMatrices({'overlap', 'coreHamilt'});
            properties.overlapMat = matrices{1};
            properties.coreHamilt = matrices{2};
            scalars = G09RSCF.G09ReadScalars({'NumElectrons', 'NucRepEnergy'});
            properties.numElectrons = scalars{1};
            properties.nucRepEnergy = scalars{2};
            properties.matpsi2 = [];
            obj@UHF(properties);
            obj.info = info;
        end
        
        function [densVec, orbital] = HarrisGuess(obj)
            info_ = obj.info;
            info_.harris = [];
            G09RSCF.RunG09(info_);
            matrices = G09RSCF.G09ReadMatrices({'HarrisGuessMOAlpha', 'HarrisGuessMOBeta'});
            orbital = matrices;
            densVec = obj.OrbToDensVec(orbital);
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(obj, orbital)
            obj.info.orbAlpha = orbital{1};
            obj.info.orbBeta = orbital{2};
            G09RSCF.RunG09(obj.info);
            matrices = G09RSCF.G09ReadMatrices({'fockAlpha', 'fockBeta'});
            fockVec = [reshape(matrices{1}, [], 1), reshape(matrices{2}, [], 1)];
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            energy = obj.SCFEnergy@G09RSCF(fockVec, densVec);
        end
        
        function energy = DampedSCFEnergy(obj, fockVec, densVec, dampingCoeff, guessOrbital)
            info_ = obj.info;
            info_.orbAlpha = guessOrbital{1};
            info_.orbBeta = guessOrbital{2};
            info_.dampingCoeff = dampingCoeff;
            G09RSCF.RunG09(info_);
            scalars = G09RSCF.G09ReadScalars({'dampedEnergy'});
            energy = scalars{1}(1);
        end
        
    end
    
end
