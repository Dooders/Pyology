import logging
from typing import TYPE_CHECKING, Tuple

from utils.command_data import CommandData
from utils.tracking import execute_command

from .common_reactions import KrebsCycleReactions
from .energy_calculations import (
    calculate_energy_state,
    calculate_total_adenine_nucleotides,
)
from .exceptions import KrebsCycleError, ReactionError
from .pathway import Pathway

if TYPE_CHECKING:
    from .organelle import Organelle

logger = logging.getLogger(__name__)


class KrebsCycle(Pathway):
    """
    Class representing the Krebs Cycle (Citric Acid Cycle) pathway.

    The Krebs Cycle is a series of chemical reactions used by all aerobic organisms
    to release stored energy through the oxidation of acetyl-CoA derived from carbohydrates,
    fats, and proteins.

    Methods
    -------
    run:
        Executes the Krebs Cycle pathway for a given number of acetyl-CoA units.
    cycle:
        Executes one complete cycle of the Krebs Cycle.
    """

    time_step = 1
    reactions = KrebsCycleReactions

    def __init__(self, debug=False):
        self.debug = debug

    def run(
        self, organelle: "Organelle", acetyl_coa_units: float, logger: logging.Logger
    ) -> Tuple[float, float, float]:
        """
        Executes the Krebs Cycle pathway for a given number of acetyl-CoA units.

        Parameters
        ----------
        organelle: Organelle
            The organelle to run the Krebs Cycle pathway on.
        acetyl_coa_units: float
            The number of acetyl-CoA units to process.
        logger: logging.Logger
            The logger to use for logging messages.

        Returns
        -------
        Tuple[float, float, float]:
            A tuple containing the final energy state, final adenine nucleotides,
            and the number of CO2 produced.

        Raises
        ------
        KrebsCycleError:
            If the Krebs Cycle pathway fails to complete.
        """
        try:
            if acetyl_coa_units <= 0:
                raise KrebsCycleError(
                    "The number of acetyl-CoA units must be positive."
                )

            initial_energy = calculate_energy_state(organelle, logger)
            initial_adenine = calculate_total_adenine_nucleotides(organelle)

            co2_produced = 0
            for i in range(int(acetyl_coa_units)):
                cycle_results = execute_command(
                    organelle,
                    CommandData(
                        obj=self,
                        command=self.cycle,
                        tracked_attributes=[
                            "ATP",
                            "ADP",
                            "AMP",
                            "Acetyl_CoA",
                            "NAD",
                            "NADH",
                            "FAD",
                            "FADH2",
                            "CO2",
                        ],
                        args=(),
                        kwargs={},
                    ),
                    logger=logger,
                )
                co2_produced += cycle_results.result

            logger.info(f"Krebs Cycle completed. Produced {co2_produced} CO2.")
            logger.info(f"Final metabolite levels: {organelle.metabolites.quantities}")

            final_energy = calculate_energy_state(organelle, logger)
            final_adenine = calculate_total_adenine_nucleotides(organelle)

            # Check energy conservation
            energy_difference = final_energy - initial_energy
            if abs(energy_difference) > 1e-6:
                logger.warning(
                    f"Energy not conserved in Krebs Cycle. Difference: {energy_difference}"
                )

            # Check adenine nucleotide conservation
            adenine_difference = final_adenine - initial_adenine
            if abs(adenine_difference) > 1e-6:
                logger.warning(
                    f"Adenine nucleotides not conserved in Krebs Cycle. Difference: {adenine_difference}"
                )

            return final_energy, final_adenine, co2_produced

        except Exception as e:
            logger.error(f"Error during Krebs Cycle: {str(e)}")
            raise KrebsCycleError(f"Krebs Cycle failed: {str(e)}")

    def cycle(self, organelle: "Organelle") -> int:
        """
        Perform one complete cycle of the Krebs Cycle.

        Parameters
        ----------
        organelle: Organelle
            The organelle to run the Krebs Cycle on.

        Returns
        -------
        int:
            The number of CO2 molecules produced in this cycle.

        Raises
        ------
        KrebsCycleError:
            If the cycle fails to complete.
        """
        logger = logging.getLogger(__name__)  # Get logger within the method
        logger.info("Starting Krebs Cycle")
        co2_produced = 0

        try:
            self.reactions.citrate_synthase.transform(organelle=organelle)
            self.reactions.aconitase.transform(organelle=organelle)
            self.reactions.isocitrate_dehydrogenase.transform(organelle=organelle)
            co2_produced += 1
            self.reactions.alpha_ketoglutarate_dehydrogenase.transform(
                organelle=organelle
            )
            co2_produced += 1
            self.reactions.succinyl_coa_synthetase.transform(organelle=organelle)
            self.reactions.succinate_dehydrogenase.transform(organelle=organelle)
            self.reactions.fumarase.transform(organelle=organelle)
            self.reactions.malate_dehydrogenase.transform(organelle=organelle)

        except ReactionError as e:
            logger.error(f"Krebs Cycle failed: {str(e)}")
            raise KrebsCycleError(f"Krebs Cycle failed: {str(e)}")

        logger.info(f"Krebs Cycle completed. CO2 produced: {co2_produced}")
        return co2_produced
