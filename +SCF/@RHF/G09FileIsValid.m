function fileIsValid = G09FileIsValid(~)
logFile = fopen('temp.log');
currLine = '';
fileIsValid = false;
while ischar(currLine)
    if ~isempty(regexp(currLine, ' Normal termination ', 'ONCE'))
        fileIsValid = true;
    end
    currLine = fgetl(logFile);
end
fclose(logFile);
end
