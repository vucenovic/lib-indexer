A=imresize(imread('alphabet\..\A.bmp'),[42,24]);B=imresize(imread('alphabet\..\B.bmp'),[42,24]);
C=imresize(imread('alphabet\..\C.bmp'),[42,24]);D=imresize(imread('alphabet\..\D.bmp'),[42,24]);
E=imresize(imread('alphabet\..\E.bmp'),[42,24]);F=imresize(imread('alphabet\..\F.bmp'),[42,24]);
G=imresize(imread('alphabet\..\G.bmp'),[42,24]);H=imresize(imread('alphabet\..\H.bmp'),[42,24]);
I=imresize(imread('alphabet\..\I.bmp'),[42,24]);J=imresize(imread('alphabet\..\J.bmp'),[42,24]);
K=imresize(imread('alphabet\..\K.bmp'),[42,24]);L=imresize(imread('alphabet\..\L.bmp'),[42,24]);
M=imresize(imread('alphabet\..\M.bmp'),[42,24]);N=imresize(imread('alphabet\..\N.bmp'),[42,24]);
O=imresize(imread('alphabet\..\O.bmp'),[42,24]);P=imresize(imread('alphabet\..\P.bmp'),[42,24]);
Q=imresize(imread('alphabet\..\Q.bmp'),[42,24]);R=imresize(imread('alphabet\..\R.bmp'),[42,24]);
S=imresize(imread('alphabet\..\S.bmp'),[42,24]);T=imresize(imread('alphabet\..\T.bmp'),[42,24]);
U=imresize(imread('alphabet\..\U.bmp'),[42,24]);V=imresize(imread('alphabet\..\V.bmp'),[42,24]);
W=imresize(imread('alphabet\..\W.bmp'),[42,24]);X=imresize(imread('alphabet\..\X.bmp'),[42,24]);
Y=imresize(imread('alphabet\..\.bmp'),[42,24]);Z=imresize(imread('alphabet\..\Z.bmp'),[42,24]);

alphabet=[A B C D E F G H I J K L M N O P Q R S T U V W X Y Z];

templates=mat2cell(alphabet,42,[24 24 24 24 24 24 24 ...
    24 24 24 24 24 24 24 ...
    24 24 24 24 24 24 24 ...
    24 24 24 24 24]);

save ('templates','templates')
