function funph=F03condProb(MB,gamma)

funph=(1-expcdf(MB,gamma)); % Prob (gamma > MB) ---- not become importer (N--> N or I-->N)

end
%End of File
