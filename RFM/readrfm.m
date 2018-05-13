function [RFMcoefficient, regularization] = readrfm( filepath )
%   readrfm 输入路径, 读取rpb文件
%   RFMcoefficient为系数矩阵，四列依次对应
%   L分子 L分母 S分子 S分母
%   L_Num L_Den S_Num S_Den
%   regularization为正则化参数
%   元素依次为
%   lineOffset = regularization(1);
%   sampOffset = regularization(2);
%   latOffset = regularization(3);
%   longOffset = regularization(4);
%   heightOffset = regularization(5);
%   lineScale = regularization(6);
%   sampScale = regularization(7);
%   latScale = regularization(8);
%   longScale = regularization(9);
%   heightScale = regularization(10);

    fid = fopen(filepath,'r');
    textscan(fid,'%*s',18);
    %正则化参数
    temp = textscan(fid,'%*s%*s%f;',10);
    regularization = temp{1};
    %L分子系数
    textscan(fid,'%*s',3);
    L_Num = fscanf(fid,'%f,',[20,1]);
    %L分母系数
    textscan(fid,'%*s',4);
    L_Den = fscanf(fid,'%f,',[20,1]);
    %S分子系数
    textscan(fid,'%*s',4);
    S_Num = fscanf(fid,'%f,',[20,1]);
    %S分母系数
    textscan(fid,'%*s',4);
    S_Den = fscanf(fid,'%f,',[20,1]);
    fclose(fid);
    RFMcoefficient = [L_Num, L_Den, S_Num, S_Den];

end

