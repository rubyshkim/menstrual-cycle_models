function y = Hill_plus(y,km,n)
	% Positive Hill Funuction

	y = ((y/km).^n)./(1+(y/km).^n);

end