imageBWD = imread('data/ZY3_TLC_E85.9_N46.1_20140804_L1A0001804812-BWD.tiff');
imageFWD = imread('data/ZY3_TLC_E85.9_N46.1_20140804_L1A0001804812-FWD.tiff');


littleImageBWD = ImageLinearFunction(imageBWD, 0.005);
littleImageFWD = ImageLinearFunction(imageFWD, 0.005);

imwrite(littleImageBWD, 'data/littleImageBWD.tiff');
imwrite(littleImageFWD, 'data/littleImageFWD.tiff');

