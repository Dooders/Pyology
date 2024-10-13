"""
Signal Transduction Pathways
    Modeling Cascades:
        Implement sequences where one activated molecule activates others.
    Feedback Loops:
        Incorporate positive and negative feedback mechanisms.
"""

from typing import Dict


class Enzyme:
    """
    Represents an enzyme that catalyzes a biochemical reaction.

    Parameters
    ----------
    name : str
        The name of the enzyme.
    vmax : float
        The maximum rate of the enzyme.
    km : float
        The Michaelis-Menten constant.
    inhibitors : Dict[str, float], optional
        Inhibitors and their inhibition constants (default is None).
    activators : Dict[str, float], optional
        Activators and their activation constants (default is None).

    Methods
    -------
    catalyze : function
        Catalyzes a biochemical reaction.
    """

    def __init__(
        self, name: str, vmax: float, km: float, inhibitors=None, activators=None
    ):
        self.name = name
        self.vmax = vmax
        self.km = km
        self.inhibitors = inhibitors if inhibitors else {}
        self.activators = activators if activators else {}

    def calculate_rate(
        self, substrate_concentration: float, metabolite_levels: Dict[str, float]
    ) -> float:
        """
        Calculate the reaction rate, considering inhibitors and activators.

        Parameters
        ----------
        substrate_concentration : float
            The concentration of the substrate.
        metabolite_levels : Dict[str, float]
            Current levels of metabolites that may inhibit or activate the enzyme.

        Returns
        -------
        float
            The reaction rate.
        """
        # Adjust km based on inhibitors
        km_effective = self.km
        for inhibitor, ki in self.inhibitors.items():
            inhibitor_concentration = metabolite_levels.get(inhibitor, 0)
            km_effective *= 1 + inhibitor_concentration / ki

        # Adjust vmax based on activators
        vmax_effective = self.vmax
        for activator, ka in self.activators.items():
            activator_concentration = metabolite_levels.get(activator, 0)
            vmax_effective *= 1 + activator_concentration / ka

        return (vmax_effective * substrate_concentration) / (
            km_effective + substrate_concentration
        )