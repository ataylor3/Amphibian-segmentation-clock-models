% USES LOCAL EXTREMA TO FIND PERIOD AND AMPLITUDE OF EXPRESSION
function osc_behavior = f(t,tfinal,protein,mRNA)

    % look for local extrema
    % islocal functions give logical array, 1s at extrema, 0s otherwise

    p_max = protein(islocalmax(protein)); % vector of local protein max
    p_min = protein(islocalmin(protein)); % vector of local protein min
    m_max = mRNA(islocalmax(mRNA)); % vector of local mRNA max
    m_min = mRNA(islocalmin(mRNA)); % vector of local mRNA min
    
    t_p_min = t(islocalmin(protein)); % vector of time stamps, local p min
    t_m_min = t(islocalmin(mRNA)); % vector of time stamps, local m min
    
    % cut out first 5 cycles of oscillation to avoid skewing data
    p_max = p_max(1,6:end);
    p_min = p_min(1,5:end);
    m_max = m_max(1,6:end);
    m_min = m_min(1,5:end);
    
    t_p_min = t_p_min(1,5:end);
    t_m_min = t_m_min(1,5:end);

    % find period of expression
    periodicity_mRNA = P(t_m_min,tfinal,m_min,m_max);
    periodicity_protein = P(t_p_min,tfinal,p_min,p_max);
    
    % find amplitude of expression
    amplitude_mRNA = Amp(m_max,m_min);
    amplitude_protein = Amp(p_max,p_min);
    
    % store results
    osc_behavior = [periodicity_mRNA, periodicity_protein, amplitude_mRNA, amplitude_protein];
end
