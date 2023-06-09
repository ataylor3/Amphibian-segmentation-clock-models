function Amp = a(max_mol,min_mol) % AMPLITUDE OF EXPRESSION

if size(max_mol,2) == size(min_mol,2) % make sure corresponding maxima and minima are matched up
    amp = max_mol - min_mol;
elseif size(max_mol,2) > size(min_mol,2)
    amp = max_mol(:,1:size(min_mol,2)) - min_mol;
elseif size(min_mol,2) > size(max_mol,2)
    amp = max_mol - min_mol(:,1:size(max_mol,2));
end


Amp = mean(amp); % take average amplitude

end

