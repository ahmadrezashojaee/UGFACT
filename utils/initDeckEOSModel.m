function model = initDeckEOSModel(deck)
% Set up a EquationOfState model from a parsed deck
%
% SYNOPSIS:
%   model = initDeckEOSModel(deck)
%
% DESCRIPTION:
%   Detailed description of function
%
% PARAMETERS:
%   deck  - Input deck with compositional keywords.
%
% RETURNS:
%   model - Initialized EquationOfState model.

    names = deck.PROPS.CNAMES;
    Tcrit = deck.PROPS.TCRIT';
    Pcrit = deck.PROPS.PCRIT';
    if isfield(deck.PROPS, 'VCRIT')
        Vcrit = deck.PROPS.VCRIT';
    else
        R = 8.3144598;
        Vcrit = R*deck.PROPS.ZCRIT'.*Tcrit./Pcrit;
    end
    acf   = deck.PROPS.ACF';
    mw    = deck.PROPS.MW';
    
    T_ref = 273.15 + 15;
    P_ref = 1*atm;
    
    fluid = CompositionalMixture(names, Tcrit, Pcrit, Vcrit, acf, mw, T_ref, P_ref);
    
    if isfield(deck.PROPS, 'BIC') && any(deck.PROPS.BIC ~= 0)
        fluid = fluid.setBinaryInteraction(deck.PROPS.BIC);
    end
    
    model = EquationOfStateModel([], fluid);
    
    eos_type = deck.PROPS.EOS;
    if iscell(eos_type) && numel(eos_type) > 1,
        warning('Multiple EOS regions not supported, using region 1');
        eos_type = eos_type{1};
    end
    if isfield(deck.PROPS, 'PRCORR') && (isempty(eos_type) || strcmp(eos_type, 'PR'))
        model = model.setType('PRCORR');
    else
        model = model.setType(eos_type);
    end
    if isfield(deck.PROPS, 'OMEGAA')
        omegaA = deck.PROPS.OMEGAA(1,:);
        defaulted = isnan(omegaA);
        switch eos_type
            case 'PR'
                val = 0.457235529;
            case {'RK', 'SRK', 'ZJ'}
                val = 0.4274802;
        end
        omegaA(defaulted) = val;
        model.omegaA = omegaA;
        % @@ sinlge value for all components for now
        model.omegaA = omegaA(1);
    end
    if isfield(deck.PROPS, 'OMEGAB')
        omegaB = deck.PROPS.OMEGAB(1,:);
        defaulted = isnan(omegaB);
        switch eos_type
            case 'PR'
                val = 0.077796074;
            case {'RK', 'SRK', 'ZJ'}
                val =  0.08664035;
        end
        omegaB(defaulted) = val;
        model.omegaB = omegaB;
        % @@ sinlge value for all components for now
        model.omegaB = omegaB(1);
    end
    if isfield(deck.PROPS, 'SSHIFT')
        model.PropertyModel.volumeShift = deck.PROPS.SSHIFT;
    end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
