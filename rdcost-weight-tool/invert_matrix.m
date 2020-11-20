A = dlmread("/dev/stdin");
B = inv(A);
dlmwrite("/dev/stdout", B, " ");
