% OPTICAL CHARACTER RECOGNITION
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% OCR algorithm, which consists of two parts:
% 1. Preprocessing
% 2. NCC or SSD template matching algorithm

function label = ocr(label)

STRATEGY = "NCC";

%patches = preprocessing(label);
templates = loadTemplates();

end

function templates = loadTemplates()
currentFolder = pwd;
pathTemplate = strcat(pwd,'\templates\');
files = dir(fullfile(pathTemplate,'*.bmp'));
templates = [];
for k = 1:numel(files)
    file = fullfile(pathTemplate, files(k).name);
    template = imread(file);
    templates = [templates, struct("template", template, "char", files(k).name)];
end
end