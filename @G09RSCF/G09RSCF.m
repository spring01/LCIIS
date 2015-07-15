classdef G09RSCF < RHF
    
    properties (Access = protected)
        
        info;
        
    end
    
    methods
        
        function obj = G09RSCF(info)
            G09RSCF.RunG09(info);
            matrices = G09RSCF.G09ReadMatrices({'overlap', 'coreHamilt'});
            properties.overlapMat = matrices{1};
            properties.coreHamilt = matrices{2};
            scalars = G09RSCF.G09ReadScalars({'NumElectrons', 'NucRepEnergy'});
            properties.numElectrons = scalars{1};
            properties.nucRepEnergy = scalars{2};
            properties.matpsi2 = [];
            obj@RHF(properties);
            obj.info = info;
        end
        
        function [densVec, orbital] = HarrisGuess(obj)
            info_ = obj.info;
            info_.harris = [];
            G09RSCF.RunG09(info_);
            matrices = G09RSCF.G09ReadMatrices({'HarrisGuessMO'});
            orbital = matrices{1};
            densVec = obj.OrbToDensVec(orbital);
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(obj, orbital)
            obj.info.orbAlpha = orbital;
            G09RSCF.RunG09(obj.info);
            matrices = G09RSCF.G09ReadMatrices({'fockAlpha'});
            fockVec = reshape(matrices{1}, [], 1);
        end
        
        function energy = SCFEnergy(obj, fockVec, densVec)
            scalars = G09RSCF.G09ReadScalars({'totalEnergy'});
            energy = scalars{1};
        end
        
        function energy = DampedSCFEnergy(obj, fockVec, densVec, dampingCoeff, guessOrbital)
            info_ = obj.info;
            info_.orbAlpha = guessOrbital;
            info_.dampingCoeff = dampingCoeff;
            G09RSCF.RunG09(info_);
            scalars = G09RSCF.G09ReadScalars({'dampedEnergy'});
            energy = scalars{1}(1);
        end
        
    end
    
    methods (Static)
        
        scalars = G09ReadScalars(types);
        matrices = G09ReadMatrices(types);
        
    end
    
    methods (Static, Access = protected)
        
        function RunG09(info)
            fileStr = G09RSCF.G09InputStr(info);
            gjfFile = fopen('temp.gjf', 'w');
            fprintf(gjfFile, '%s', fileStr);
            fclose(gjfFile);
            system('g09 temp.gjf');
            if(~G09RSCF.G09FileIsValid())
                throw(MException('G09RSCF:RunG09', 'g09 did not terminate correctly'));
            end
        end
        
    end
    
    methods (Static, Access = private)
        
        fileStr = G09InputStr(info);
        fileIsValid = G09FileIsValid();
        
    end
    
end
