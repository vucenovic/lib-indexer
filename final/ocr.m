% OPTICAL CHARACTER RECOGNITION
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% OCR algorithm, which consists of two parts:
% 1. Preprocessing
% 2. NCC or SSD template matching algorithm

function label = ocr(label)

STRATEGY = "NCC";
label = imread('label_1.png');
patches = preprocessing(label);
templates = loadTemplates();

for k = 1:numel(patches)
    values = []
    wordOne = '';
    wordTwo = '';
    wordThree = '';
    for j = 1:numel(templates)
        if STRATEGY == "NCC";
            value = ncc(templates(j).template, patches(k).image);
            values = [values, value];
            [bestCorr, index] = max(values, [], 'all');
            matchedChar = templates(index).char;
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