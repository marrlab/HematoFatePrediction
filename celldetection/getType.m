function [type,onset,typeString] = getType(celltrack)
    % wavelength1: PU.1
    % wavelength2: GATA (MEP) -> type 1
    % wavelength3: FCgamma(GMP) -> type 2
    % nothing -> type 0
    if ~isempty(find(celltrack.wavelength_2, 1)) && isempty(find(celltrack.wavelength_3, 1))
        type = 1;
        onset = find(celltrack.wavelength_2,1,'first');
        typeString = 'MEP';
    elseif ~isempty(find(celltrack.wavelength_3, 1)) && isempty(find(celltrack.wavelength_2, 1))
        type = 2;
        onset = find(celltrack.wavelength_3,1,'first');
        typeString = 'GMP';
    elseif ~isempty(find(celltrack.wavelength_2, 1)) && ~isempty(find(celltrack.wavelength_3, 1))
        type = -1;
        onset = -1;
        typeString = 'BOTH';
    else
        type = 0;
        onset = -1;
        typeString = 'UNKNOWN';
    end

end