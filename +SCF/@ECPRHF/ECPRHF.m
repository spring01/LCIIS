classdef ECPRHF < SCF.RHF
    
    methods
        
        function obj = ECPRHF(properties)
            obj@SCF.RHF(properties);
        end
        
    end
    
    methods (Static)
        
        prop = InfoInterface(info);
        
    end
    
    methods (Static, Access = protected)
        
        ecpMat = G09ReadECPMatrix();
        
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
    
    methods (Static, Access = private)
        
        fileStr = G09InputStr(info);
        fileIsValid = G09FileIsValid();
        
    end
    
end
