function y = Hill_minus(y,km,n)
	%Negative Hill Function

	y = 1./(1+(y/km).^n);

end