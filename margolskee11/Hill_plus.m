function y = Hill_plus(y,km,n)
	% Positive Hill Function

	y = ((y/km).^n)./(1+(y/km).^n);

end