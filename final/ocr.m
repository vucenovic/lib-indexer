% OPTICAL CHARACTER RECOGNITION
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% OCR algorithm, which consists of two parts:
% 1. Preprocessing
% 2. NCC or SSD template matching algorithm

function label = ocr(label)

STRATEGY = "SSD";
label = imread('label_1.png');
patch = preprocessing(label);
templates = loadTemplates();

for k = 1:numel(patch)
    wordOne = '';
    wordTwo = '';
    wordThree = '';
    bestCorr = 10.0;
    for j = 1:numel(templates)
        if STRATEGY == "NCC";
            corrMatrix = ncc(templates(j).template, patch(k).image);
            value = max(corrMatrix, [], 'all');
            if value > bestCorr
                bestCorr = value;
                bestIndex = j;
            matchedChar = templates(bestIndex).char;
            end
        elseif STRATEGY == "SSD";
            value = ssd_naive(templates(j).template, patch(k).image);
            if value < bestCorr
                bestCorr = value;
                bestIndex = j;
            matchedChar = templates(bestIndex).char;
            end
        end
    end
end


end

function templates = loadTemplates()
currentFolder = pwd;
pathTemplate = strcat(pwd,'\templates\');
files = dir(fullfile(pathTemplate,'*.bmp'));
templates = [];
for k = 1:numel(files)
    file = fullfile(pathTemplate, files(k).name);
    template = imread(file);
    charname = extractBefore(files(k).name, ".");
    templates = [templates, struct("template", template, "char", charname)];
end
end