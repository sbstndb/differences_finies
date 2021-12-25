SA = sprand(8, 8, 0.2);
SA = full(SA);
[A, CA, RA] = myldlt_to_delete(SA);
[sss1, sss2] = size(SA)
x = ones(sss2, 1);

y = mydspmv(A, CA, RA, x');