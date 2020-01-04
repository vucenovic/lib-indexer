% OPTICAL CHARACTER RECOGNITION
% AUTHOR: ALEKSANDAR VUCENOVIC, 01635282

% OCR algorithm, which consists of two parts:
% 1. Preprocessing
% 2. NCC or SSD template matching algorithm

function label = ocr(label)

STRATEGY = "NCC";

patches = preprocessing(label);
for j = 1:length(patches)
    if(STRATEGY=="NCC")
       ncc(); 
    elseif(STRATEGY=="SSD")
       
    end
end





end