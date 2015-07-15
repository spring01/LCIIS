function scalars = G09ReadScalars(types)
keyword = cell(1, length(types));
for iSca = 1:length(types)
    if(strcmpi(types{iSca}, 'numElectrons'))
        keyword{iSca} = ' alpha electrons ';
    elseif(strcmpi(types{iSca}, 'totalEnergy'))
        keyword{iSca} = ' E= ';
    elseif(strcmpi(types{iSca}, 'nucRepEnergy'))
        keyword{iSca} = ' nuclear repulsion energy ';
    elseif(strcmpi(types{iSca}, 'dampedEnergy'))
        keyword{iSca} = ' Delta-E= ';
    end
end

scalars = cell(1, length(types));
logFile = fopen('temp.log');
currLine = '';
while(ischar(currLine))
    for iSca = 1:length(types)
        if(~isempty(regexp(currLine, keyword{iSca}, 'ONCE')))
            scaCell = regexp(currLine, '([+-]?[0-9]+(\.[0-9]+)?(D[+-][0-9]+)?)', 'match');
            scalars{iSca} = str2num(char(scaCell{:})); %#ok
            break;
        end
    end
    currLine = fgetl(logFile);
end
fclose(logFile);
end
