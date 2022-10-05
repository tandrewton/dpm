function data = readDataBlocks(fstr, numCols)

fileID = fopen(fstr);
block=1;

while (~feof(fileID))
    %fprintf('block: %s\n', num2str(block))
    formatString = repmat('%f',1, numCols);
    inputText = textscan(fileID, formatString);
    data{block,1} = cell2mat(inputText);
    [numRows, numCols] = size(data{block});
    eob = textscan(fileID, '%s', 1, 'delimiter', '\n');
    block = block+1;
end