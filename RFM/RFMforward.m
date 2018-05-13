function SL = RFMforward(BLH, RFMcoefficient, regularization)
%   RFMforward 利用RFM模型正算
%   输入经纬度高程以, 正则化向量, RFM模型参数计算像点坐标
%   函数内部进行正则化，即输入Lat Lon Height(不是V U W)
%   直接得到像点坐标(S, L)，即在外部不用再正则化以及反解
%   输出: S = SL(1) L = SL(2)
%   输入: Lat = BLH(1) Lon = BLH(2) Height = BLH(3)
%   RFMcoefficient与regularization的值直接输函数readrfm返回值即可
    lineOffset = regularization(1);
    sampOffset = regularization(2);
    latOffset = regularization(3);
    longOffset = regularization(4);
    heightOffset = regularization(5);
    lineScale = regularization(6);
    sampScale = regularization(7);
    latScale = regularization(8);
    longScale = regularization(9);
    heightScale = regularization(10);

    Lat = BLH(1); 
    Lon = BLH(2); 
    Height = BLH(3);
    
    L_Num = RFMcoefficient(:,1);
    L_Den = RFMcoefficient(:,2);
    S_Num = RFMcoefficient(:,3);
    S_Den = RFMcoefficient(:,4);

    %正则化BLH
    P = (Lat - latOffset) / latScale;%V
    L = (Lon - longOffset) / longScale;%U
    H = (Height - heightOffset) / heightScale;%W
    
    %rfm计算
    X = ones(1,20);
    X(2) = L; X(3) = P; X(4) = H;
    X(5) = L * P; X(6) = L * H; X(7) = P * H;
    X(8) = L * L; X(9) = P * P; X(10) = H * H;
    X(11) = P * L * H; X(12) = L * L * L; X(13) = L * P * P;
    X(14) = L * H * H; X(15) = L * L * P; X(16) = P * P * P;
    X(17) = P * H * H; X(18) = L * L * H; X(19) = P * P * H;
    X(20) = H * H * H;
    
    x = (X * S_Num) / (X * S_Den);
    y = (X * L_Num) / (X * L_Den);
    
    %反正则化
    SL1(1) = x * sampScale + sampOffset;
    SL1(2) = y * lineScale + lineOffset;
    
    if(length(regularization)==16)
       SL(1) = SL1(1) + (regularization(14)+regularization(15)*SL1(1) + regularization(16)*SL1(2));
       SL(2) = SL1(2) + (regularization(11)+regularization(12)*SL1(1) + regularization(13)*SL1(2));
    else
        SL(1) = SL1(1);
        SL(2) = SL1(2);
    end
    
end

