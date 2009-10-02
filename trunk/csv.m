function csv(A, format = "%f")
% csv(A, format) prints a matrix to stdout as formatted, comma-separated values.
% The format must be suitable for passing to printf and defaults to "%f".
for i = 1:size(A)(1)
	for j = 1:size(A)(2)
		if j < size(A)(2)
			f = cstrcat(format, ", ");
		else
			if i < size(A)(1)
				f = cstrcat(format, ",");
			else
				f = format;
			endif
		endif
		printf(f, A(i, j))
	endfor
	printf("\n")
endfor
