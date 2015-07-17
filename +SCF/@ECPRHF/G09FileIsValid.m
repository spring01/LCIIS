function fileIsValid = G09FileIsValid()
logFile = fopen('temp.log');
currLine = '';
fileIsValid = 0;
while(ischar(currLine))
    if(~isempty(regexp(currLine, ' Normal termination ', 'ONCE')))
        fileIsValid = 1;
    end
    currLine = fgetl(logFile);
end
fclose(logFile);
end
