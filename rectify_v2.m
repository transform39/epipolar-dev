% 读取RFM参数
[RFMcoefficientBWD, regularizationBWD] = readrfm('ZY3_TLC_E85.9_N46.1_20140804_L1A0001804812-BWD.rpb'); % 左片
[RFMcoefficientFWD, regularizationFWD] = readrfm('ZY3_TLC_E85.9_N46.1_20140804_L1A0001804812-FWD.rpb'); % 右片

% 测区平均航高
avgHeight = 505984; 
deltaHeight = 10;
upHeight = avgHeight + deltaHeight;
downHeight = avgHeight - deltaHeight;

% 影像读取
% imageBWD = imread('data/ZY3_TLC_E85.9_N46.1_20140804_L1A0001804812-BWD.tiff');
% imageFWD = imread('data/ZY3_TLC_E85.9_N46.1_20140804_L1A0001804812-FWD.tiff');


% littleImageBWD = ImageLinearFunction(imageBWD, 0.005);
% littleImageFWD = ImageLinearFunction(imageFWD, 0.005);

littleImageBWD = imread('littleImageBWD.tiff');
littleImageFWD = imread('littleImageFWD.tiff');


% 影像尺寸
[rowsBWD, colsBWD] = size(littleImageBWD);
[rowsFWD, colsFWD] = size(littleImageFWD);

imageEpi1 = uint8(zeros(1.5*rowsBWD, 1.5*colsBWD));
imageEpi2 = uint8(zeros(1.5*rowsBWD, 1.5*colsBWD));


%测试方程系数集
ParameterAB1 = zeros(rowsBWD + colsBWD - 1, 2);
ParameterAB2 = zeros(rowsFWD + colsFWD - 1, 2);


for m = 1:colsBWD
    
    % ----------步骤1：左片正算计算点1和点2的大地坐标----------
    Sa = 1;
    La = m;
    
    % 组成点1和点2的SLH量，便于调用RFMreverse函数
    SLH1 = [Sa, La, upHeight];
    SLH2 = [Sa, La, downHeight];

    % 计算点1和点2的大地坐标，都用左片正算
    BL1 = RFMreverse(SLH1, RFMcoefficientBWD, regularizationBWD);
    BL2 = RFMreverse(SLH2, RFMcoefficientBWD, regularizationBWD);

    % ----------步骤2：右片反算计算点b和点c的影像坐标----------

    % 组成点1和点2的BLH量，便于调用RFMforward函数
    BLH1 = [BL1, upHeight];
    BLH2 = [BL2, downHeight];

    % 计算点b和点c的影像坐标，都用右片算
    SLb = RFMforward(BLH1, RFMcoefficientFWD, regularizationFWD);
    SLc = RFMforward(BLH2, RFMcoefficientFWD, regularizationFWD);
    
    Sb = SLb(1);
    Sc = SLc(1);
    Lb = SLb(2);
    Lc = SLc(2);

    % ----------步骤3：右片正算计算点3的大地坐标----------

    SLH3 = [SLc, upHeight];

    BL3 = RFMreverse(SLH3, RFMcoefficientFWD, regularizationFWD);

    % ----------步骤4：左片反算计算点d的影像坐标----------

    BLH3 = [BL3, upHeight];

    SLd = RFMforward(BLH3, RFMcoefficientBWD, regularizationBWD);

    Sd = SLd(1);
    Ld = SLd(2);

    % ----------步骤5：计算直线方程----------

    if (La -Ld ~= 0) 
        A1 = (Sa - Sd) / (La - Ld);
    end
    B1 = Sa - A1 * La;
    if (Lb -Lc ~= 0) 
        A2 = (Sb - Sc) / (Lb - Lc);
    end
    B2 = Sb - A2 * Lb;
    
    ParameterAB1(m, 1) = A1;
    ParameterAB1(m, 2) = B1;
    ParameterAB2(m, 1) = A2;
    ParameterAB2(m, 2) = B2;
 
    deltaBY = 1 / sqrt(1 + (A1 * A1));
    
    % 每行
    indexBWD = 1;
    for stepY = 1:deltaBY:rowsBWD - 1
        stepX = A1*stepY - B1;
        if (stepX > colsBWD - 1 || floor(stepY) <= 0 || floor(stepX) <=0)
           continue 
        end
        fuckX = stepX - floor(stepX);
        fuckY = stepY - floor(stepY);
        value = (1 - fuckX)*(1- fuckY)*littleImageBWD(floor(stepY), floor(stepX)) +...
        (1 - fuckX)*fuckY*littleImageBWD(floor(stepY), floor(stepX)+1) + ...
        fuckX*(1 - fuckY)*littleImageBWD(floor(stepY) + 1, floor(stepX)) + ...
        fuckX*fuckY*littleImageBWD(floor(stepY) + 1, floor(stepX) + 1);
        
        imageEpi1(m, indexBWD) = value;
        
        indexBWD = indexBWD + 1;     
    end
    
    
    deltaFY = 1 / sqrt(1 + (A2 * A2));
    
    indexFWD = 1;
    for stepY = 1:deltaFY:rowsFWD - 1
        stepX = A2 * stepY + B2;
        if (stepX > colsFWD - 1 || floor(stepY) <= 0 || floor(stepX) <=0)
           continue 
        end
        fuckX = stepX - floor(stepX);
        fuckY = stepY - floor(stepY);
        value = (1 - fuckX)*(1- fuckY)*littleImageFWD(floor(stepY), floor(stepX)) +...
        (1 - fuckX)*fuckY*littleImageFWD(floor(stepY), floor(stepX)+1) + ...
        fuckX*(1 - fuckY)*littleImageFWD(floor(stepY) + 1, floor(stepX)) + ...
        fuckX*fuckY*littleImageFWD(floor(stepY) + 1, floor(stepX) + 1);
        
        imageEpi2(m, indexFWD) = value;
        
        indexFWD = indexFWD + 1;
           
    end
    
end



