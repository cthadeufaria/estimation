

input_data = load("/home/carlos/Documents/FEUP/Semestre 3 - 2023-2024/ESD/Conteudos/Lab Activities/Lab Activity2 Data/actSheet2Data.mat");

n0 = 0;
n1 = 1;
n2 = 2;
n3 = 3;
n4 = 4;
n5 = 5;

[Y0, Phi0, sys0] = idFIR(input_data, n0);
[Y1, Phi1, sys1] = idFIR(input_data, n1);
[Y2, Phi2, sys2] = idFIR(input_data, n2);
[Y3, Phi3, sys3] = idFIR(input_data, n3);
[Y4, Phi4, sys4] = idFIR(input_data, n4);
[Y5, Phi5, sys5] = idFIR(input_data, n5);


