
function P = p(min_time_vector,tfinal,min_mol,max_mol) % PERIOD OF EXPRESSION

    P = mean(diff(min_time_vector)); % period of expression is the average number of minutes 
    % between successive local minima 

    % if the amplitude of expression < 10 for any cycle of expression, set
    % period to infinity (i.e. damped oscillations)

    % set number of complete cycles of oscillation to s
    if size(max_mol,2) == size(min_mol,2) | size(min_mol,2) > size(max_mol,2)
        s = size(max_mol,2);
    elseif size(max_mol,2) > size(min_mol,2)
        s = size(min_mol,2);
    end
    
    % go through each cycle of oscillation to check for amplitude < 10 
    for i = 1:s
        if abs(max_mol(1,i) - min_mol(1,i)) < 10
            P = inf;
            break
        end
    end
end