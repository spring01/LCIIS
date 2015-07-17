classdef ECPRHF < SCF.RHF
    
    methods
        
        function obj = ECPRHF(properties)
            obj@SCF.RHF(properties);
        end
        
        function [densVec, orbital] = HarrisGuess(obj)
            orbital = SCF.ECPRHF.G09ReadMatrix('HarrisGuessMOAlpha');
            order = SCF.ECPRHF.G09ToPsi4BasisOrder(obj.matpsi2.BasisSet_ShellNumFunctions());
            orbital = orbital(order, :);
            densVec = obj.OrbToDensVec(orbital);
        end
        
    end
    
    methods (Static)
        
        prop = InfoInterface(info);
        
    end
    
    methods (Static, Access = protected)
        
        matrix = G09ReadMatrix(type);
        
        function RunG09(info)
            fileStr = SCF.ECPRHF.G09InputStr(info);
            gjfFile = fopen('temp.gjf', 'w');
            fprintf(gjfFile, '%s', fileStr);
            fclose(gjfFile);
            system('g09 temp.gjf');
            if(~SCF.ECPRHF.G09FileIsValid())
                throw(MException('ECPRHF:RunG09', 'g09 did not terminate correctly'));
            end
        end
        
    end
    
    methods (Static, Access = protected)
        
        function order = G09ToPsi4BasisOrder(shellNfuncs)
            shell2startFunc = cumsum([1 shellNfuncs]);
            shell2startFunc = shell2startFunc(1:end-1);
            order = 1:sum(shellNfuncs);
            for i = 1:length(shellNfuncs)
                if(shellNfuncs(i) == 3) % 3p; need to change [z x y] -> [x y z]
                    order((1:3)+shell2startFunc(i)-1) = ...
                        order([3 1 2]+shell2startFunc(i)-1);
                end
            end
        end
        
    end
    
    methods (Static, Access = private)
        
        fileStr = G09InputStr(info);
        fileIsValid = G09FileIsValid();
        
    end
    
end
