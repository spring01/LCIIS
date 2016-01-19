classdef RHF < handle
    
    properties (Access = protected)
        
        overlapMat;
        invSqrtmS;
        coreHamilt;
        nucRepEnergy;
        numElectrons;
        matpsi;
        
        maxSCFIter = 200;
        RMSDensityThreshold = 1e-8;
        MaxDensityThreshold = 1e-6;
        EnergyThreshold = 1e-6;
        
    end
    
    properties (SetAccess = private)
        
        info;
        
    end
    
    methods
        
        function self = RHF(info)
            % info.charge
            % info.multiplicity
            % info.cartesian (with atomic numbers and coordinates in Angstrom)
            % info.basisSet
            % info.g09BasisSetKeyword (optional)
            % info.ecpFile (optional)
            self.info = info;
            self.matpsi = MatPsi2(info.cartesian, info.basisSet, ...
                                  info.charge, info.multiplicity);
            self.overlapMat = self.matpsi.Integrals_Overlap();
            self.invSqrtmS = inv(sqrtm(self.overlapMat));
            self.coreHamilt = self.matpsi.Integrals_Kinetic() ...
                            + self.matpsi.Integrals_Potential();
            self.nucRepEnergy = self.matpsi.Molecule_NucRepEnergy();
            totalNumElectrons = self.matpsi.Molecule_NumElectrons();
            self.numElectrons = SCF.RHF.CalcNumElectrons(totalNumElectrons, ...
                                                         info.multiplicity);
        end
        
        function RunG09(self)
            fileStr = self.G09InputStr();
            gjfFile = fopen('temp.gjf', 'w');
            fprintf(gjfFile, '%s', fileStr);
            fclose(gjfFile);
            system('g09 temp.gjf');
            if ~self.G09FileIsValid()
                throw(MException('RHF:RunG09', 'g09 failed to run'));
            end
        end
        
        function [densVec, orbital] = CoreGuess(self)
            coreFockVec = reshape(self.coreHamilt, [], 1);
            orbital = self.SolveFockVec(coreFockVec);
            densVec = self.OrbToDensVec(orbital);
        end
        
        function [densVec, orbital] = HarrisGuess(self)
            orbital = self.G09ReadMatrix('HarrisGuessMOAlpha');
            orbital = self.G09OrbToPsi4Orb(orbital);
            densVec = self.OrbToDensVec(orbital);
        end
        
    end
    
    methods (Access = protected)
        
        function densVec = OrbToDensVec(self, orbital)
            occOrb = orbital(:, 1:self.numElectrons(1));
            densVec = reshape(occOrb * occOrb', [], 1);
        end
        
        function fockVec = OrbToFockVec(self, orbital)
            occOrb = orbital(:, 1:self.numElectrons(1));
            self.matpsi.JK_CalcAllFromOccOrb(occOrb);
            gMat = 2 .* self.matpsi.JK_RetrieveJ() - self.matpsi.JK_RetrieveK();
            fockVec = reshape(self.coreHamilt + gMat, [], 1);
        end
        
        function fockVec = DensVecToFockVec(self, densVec)
            densMat = reshape(densVec, size(self.overlapMat));
            self.matpsi.JK_CalcAllFromDens(densMat);
            gMat = 2 .* self.matpsi.JK_RetrieveJ() - self.matpsi.JK_RetrieveK();
            fockVec = reshape(self.coreHamilt + gMat, [], 1);
        end
        
        function [orbital, orbEigValues] = SolveFockVec(self, fockVec)
            fockMat = reshape(fockVec, size(self.overlapMat));
            orthoFockMat = self.invSqrtmS * fockMat * self.invSqrtmS;
            orthoFockMat = (orthoFockMat + orthoFockMat') ./ 2;
            [orbitalOtho, orbEigValues] = eig(orthoFockMat);
            [orbEigValues, ascendOrder] = sort(diag(orbEigValues));
            orbitalOtho = orbitalOtho(:, ascendOrder);
            orbital = self.invSqrtmS * orbitalOtho;
        end
        
        function energy = SCFEnergy(self, fockVec, densVec)
            energy = (reshape(self.coreHamilt, [], 1) + fockVec)' * densVec ...
                   + self.nucRepEnergy;
        end
        
        function cdiis = CDIIS(self, numVectors)
            cdiis = SCF.CDIIS(self.overlapMat, numVectors, 'r');
        end
        
        function ediis = EDIIS(self, numVectors)
            ediis = SCF.EDIIS(self.coreHamilt, numVectors, 'r');
        end
        
        function adiis = ADIIS(self, numVectors)
            adiis = SCF.ADIIS(self.coreHamilt, numVectors, 'r');
        end
        
        function lciis = LCIIS(self, numVectors)
            lciis = SCF.RLCIIS(self.overlapMat, numVectors);
        end
        
        orb = G09OrbToPsi4Orb(self, orb);
        
    end
    
    methods (Access = private)
        
        fileStr = G09InputStr(self);
        fileIsValid = G09FileIsValid(self);
        matrix = G09ReadMatrix(self, type)
        
    end
    
    methods(Static, Access = private)
        
        function numElectronsAB = CalcNumElectrons(numElectronsTotal, mult)
            numElectronsA = (numElectronsTotal + mult - 1) / 2;
            if rem(numElectronsA, 1) ~= 0
                throw(MException('RHF:CalcNumElectrons', ...
                    'Number of electrons and multiplicity do not agree'));
            end
            numElectronsB = numElectronsTotal - numElectronsA;
            numElectronsAB = [numElectronsA; numElectronsB];
        end
        
    end
    
end