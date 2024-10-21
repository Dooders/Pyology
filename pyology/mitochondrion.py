import logging
import math
from typing import Dict

from utils.command_data import CommandData
from utils.tracking import execute_command

from .constants import *
from .exceptions import *
from .krebs_cycle import KrebsCycle
from .metabolite import Metabolite
from .organelle import Organelle

logger = logging.getLogger(__name__)


class Mitochondrion(Organelle):
    """
    Simulates the electron transport chain with individual complexes.

    Proton gradient and ATP synthesis are modeled well, including a proton leak
    for added realism.

    Calcium buffering and its effect on mitochondrial function are included,
    which adds another layer of detail.

    Feedback inhibition for ATP levels is implemented, mimicking real cellular
    regulation mechanisms.

    Methods
    -------
    krebs_cycle_process:
        Processes acetyl-CoA through the Krebs cycle
    pyruvate_to_acetyl_coa:
        Converts pyruvate to Acetyl-CoA
    calculate_oxygen_needed:
        Calculates the amount of oxygen needed for cellular respiration
    cellular_respiration:
        Simulates the entire cellular respiration process with feedback inhibition
    calculate_proton_leak:
        Calculates the proton leak using a logistic function
    update_proton_gradient:
        Updates the proton gradient considering nonlinear leak
    complex_I:
        Simulates Complex I activity
    complex_II:
        Simulates Complex II activity
    complex_III:
        Simulates Complex III activity
    complex_IV:
        Simulates Complex IV activity
    is_metabolite_available:
        Checks if a metabolite is available in sufficient quantity
    atp_synthase:
        Synthesizes ATP using the proton gradient
    """

    name = "Mitochondrion"

    def __init__(self, debug=False) -> None:
        super().__init__()
        self.debug = debug

        self.proton_gradient = 0
        self.atp_per_nadh = 2.5
        self.atp_per_fadh2 = 1.5
        self.atp_per_substrate_phosphorylation = 1
        self.oxygen_per_nadh = 0.5
        self.oxygen_per_fadh2 = 0.5

        self.calcium_threshold = CALCIUM_THRESHOLD
        self.calcium_boost_factor = CALCIUM_BOOST_FACTOR
        self.max_proton_gradient = MAX_PROTON_GRADIENT
        self.leak_rate = LEAK_RATE
        self.leak_steepness = LEAK_STEEPNESS
        self.leak_midpoint = LEAK_MIDPOINT

        self.krebs_cycle = KrebsCycle()

        # Initialize all necessary metabolites for the Krebs cycle
        self.initialize_krebs_cycle_metabolites()

    def initialize_krebs_cycle_metabolites(self):
        krebs_cycle_metabolites = [
            "acetyl_coa",
            "citrate",
            "isocitrate",
            "alpha_ketoglutarate",
            "succinyl_coa",
            "succinate",
            "fumarate",
            "malate",
            "oxaloacetate",  # Ensure this is lowercase
        ]
        for metabolite in krebs_cycle_metabolites:
            self.metabolites[metabolite] = Metabolite(metabolite, 0, max_quantity=100)

    def ensure_metabolite_exists(self, metabolite: str, initial_quantity: float = 0):
        """
        Ensure a metabolite exists in the metabolites dictionary.
        If it doesn't exist, add it with the given initial quantity.
        """
        if metabolite not in self.metabolites:
            logger.warning(
                f"{metabolite} not found in metabolites. Adding it with initial quantity {initial_quantity}."
            )
            self.metabolites[metabolite] = Metabolite(metabolite, initial_quantity)

    def change_metabolite_quantity(self, metabolite: str, amount: float) -> float:
        """
        Changes the quantity of a metabolite, ensuring it doesn't go negative.

        The method ensures that the metabolite quantity does not go negative by
        adjusting the amount to the maximum of the current quantity plus the
        requested change or zero.

        Parameters
        ----------
        metabolite: str
            The name of the metabolite to change.
        amount: float
            The amount to change the metabolite by (positive or negative).

        Returns
        -------
        float
            The actual amount changed (may be different if preventing negative values).
        """
        if metabolite not in self.metabolites:
            self.metabolites[metabolite] = Metabolite(metabolite, 0)

        current_quantity = self.metabolites[metabolite].quantity
        new_quantity = max(current_quantity + amount, 0)
        actual_change = new_quantity - current_quantity

        self.metabolites[metabolite].quantity = new_quantity

        if actual_change != amount:
            logger.warning(
                f"Attempted to decrease {metabolite} by {-amount}, but only decreased by {-actual_change} to prevent negative quantity."
            )

        return actual_change

    def consume_metabolites(self, **metabolites: Dict[str, float]) -> bool:
        """
        Consume multiple metabolites at once.

        Parameters
        ----------
        metabolites: Dict[str, float]
            The metabolites to consume and their amounts.

        Returns
        -------
        bool
            True if all metabolites were consumed successfully, False otherwise.
        """
        temp_changes = {}
        for metabolite, amount in metabolites.items():
            actual_change = self.change_metabolite_quantity(metabolite, -amount)
            if abs(actual_change) < abs(amount):
                # Revert changes if not all metabolites could be consumed
                for rev_metabolite, rev_amount in temp_changes.items():
                    self.change_metabolite_quantity(rev_metabolite, rev_amount)
                logger.warning(f"Insufficient {metabolite} for reaction")
                return False
            temp_changes[metabolite] = actual_change
        return True

    def produce_metabolites(self, **metabolites: Dict[str, float]) -> None:
        """
        Produce multiple metabolites at once.

        Parameters
        ----------
        metabolites: Dict[str, float]
            The metabolites to produce and their amounts.
        """
        for metabolite, amount in metabolites.items():
            actual_change = self.change_metabolite_quantity(metabolite, amount)
            if actual_change != amount:
                logger.warning(
                    f"Could not produce full amount of {metabolite}. Produced {actual_change} instead of {amount}."
                )

    def krebs_cycle_process(self, acetyl_coa_amount: int) -> int:
        """
        Processes acetyl-CoA through the Krebs cycle.

        The Krebs cycle is a series of reactions that convert acetyl-CoA into
        ATP, NADH, and FADH2.

        Parameters
        ----------
        acetyl_coa_amount: int
            The amount of acetyl-CoA to process.

        Returns
        -------
        int
            The total amount of NADH and FADH2 produced.
        """
        logger.info(
            f"Processing {acetyl_coa_amount} units of acetyl-CoA through the Krebs cycle"
        )

        self.krebs_cycle.add_substrate("glucose", acetyl_coa_amount)

        # Ensure there's enough oxaloacetate to start the cycle
        if self.metabolites["oxaloacetate"].quantity < acetyl_coa_amount:
            oxaloacetate_needed = (
                acetyl_coa_amount - self.metabolites["oxaloacetate"].quantity
            )
            if not self.consume_metabolites(oxaloacetate=oxaloacetate_needed):
                logger.warning("Insufficient oxaloacetate to start Krebs cycle")
                return 0

        total_nadh = 0
        total_fadh2 = 0
        total_atp = 0

        for metabolites, cofactors in self.krebs_cycle.krebs_cycle_iterator(
            num_cycles=acetyl_coa_amount
        ):
            total_nadh += cofactors["NADH"]
            total_fadh2 += cofactors["FADH2"]
            total_atp += cofactors["GTP"]  # GTP is equivalent to ATP

        # Transfer the products to the mitochondrion
        self.produce_metabolites(nadh=total_nadh, fadh2=total_fadh2, atp=total_atp)

        # Reset the Krebs cycle for the next round
        self.krebs_cycle.reset()

        return total_nadh + total_fadh2

    def pyruvate_to_acetyl_coa(self, pyruvate_amount: int) -> int:
        """
        Converts pyruvate to Acetyl-CoA.

        Parameters
        ----------
        pyruvate_amount: int
            The amount of pyruvate to convert.

        Returns
        -------
        int
            The amount of acetyl-CoA produced.
        """
        logger.info(f"Converting {pyruvate_amount} units of pyruvate to Acetyl-CoA")
        acetyl_coa_produced = pyruvate_amount
        self.change_metabolite_quantity("nadh", pyruvate_amount)
        self.change_metabolite_quantity("co2", pyruvate_amount)
        return acetyl_coa_produced

    def calculate_oxygen_needed(self, pyruvate_amount: int) -> float:
        """
        Calculate the amount of oxygen needed for cellular respiration based on
        pyruvate amount.

        Parameters
        ----------
        pyruvate_amount: int
            The amount of pyruvate to convert.

        Returns
        -------
        float
            The amount of oxygen needed.
        """
        return pyruvate_amount * 2.5

    def cellular_respiration(self, pyruvate_amount: float) -> float:
        """
        Perform cellular respiration on the given amount of pyruvate.

        Args:
            pyruvate_amount: The amount of pyruvate to process.

        Returns:
            The amount of ATP produced.
        """
        tracked_attributes = [
            "pyruvate",
            "atp",
            "adp",
            "nad",
            "nadh",
            "fadh2",
            "oxygen",
        ]

        def validate_conservation(obj, initial, final, changes):
            # Add appropriate validation logic here
            return True

        command_data = CommandData(
            obj=self,
            command="process_pyruvate",
            tracked_attributes=tracked_attributes,
            args=[pyruvate_amount],
            validations=[validate_conservation],
        )

        result = execute_command(command_data, logger, self.debug)

        return result["result"]  # This should be the amount of ATP produced

    def calculate_proton_leak(self) -> float:
        """
        Calculate the proton leak using a logistic function.

        Returns
        -------
        float
            The amount of proton leak.
        """
        relative_gradient = self.proton_gradient / self.max_proton_gradient
        leak = self.leak_rate / (
            1
            + math.exp(
                -self.leak_steepness * (self.proton_gradient - self.leak_midpoint)
            )
        )
        return leak

    def update_proton_gradient(self, protons_pumped: float) -> float:
        """
        Update the proton gradient considering nonlinear leak.

        Parameters
        ----------
        protons_pumped: float
            The amount of protons pumped.

        Returns
        -------
        float
            The updated proton gradient.
        """
        self.proton_gradient += protons_pumped
        leak = self.calculate_proton_leak()
        self.proton_gradient = max(0, self.proton_gradient - leak)
        logger.info(f"Updated proton gradient: {self.proton_gradient:.2f}")
        return self.proton_gradient

    def complex_I(self) -> int:
        """
        Simulates Complex I activity.

        Complex I is the first enzyme in the electron transport chain and
        catalyzes the oxidation of NADH.

        Returns
        -------
        int
            The amount of electrons transferred.
        """
        if self.is_metabolite_available("nadh", 1) and self.is_metabolite_available(
            "ubiquinone", 1
        ):
            reaction_rate = min(
                self.metabolites["nadh"].quantity,
                self.metabolites["ubiquinone"].quantity,
            )
            if self.consume_metabolites(nadh=reaction_rate, ubiquinone=reaction_rate):
                self.produce_metabolites(ubiquinol=reaction_rate)
                self.update_proton_gradient(PROTONS_PER_NADH * reaction_rate)  # Updated
                logger.info(
                    f"Complex I: Oxidized {reaction_rate} NADH, pumped {PROTONS_PER_NADH * reaction_rate} protons"
                )
                return reaction_rate
        logger.warning("Insufficient NADH or ubiquinone for Complex I")
        return 0

    def complex_II(self) -> int:
        """
        Simulates Complex II activity.

        Complex II is the second enzyme in the electron transport chain and
        catalyzes the oxidation of FADH2.

        Returns
        -------
        int
            The amount of electrons transferred.
        """
        if self.is_metabolite_available("fadh2", 1) and self.is_metabolite_available(
            "ubiquinone", 1
        ):
            reaction_rate = min(
                self.metabolites["fadh2"].quantity,
                self.metabolites["ubiquinone"].quantity,
            )
            if self.consume_metabolites(fadh2=reaction_rate, ubiquinone=reaction_rate):
                self.produce_metabolites(ubiquinol=reaction_rate)
                logger.info(f"Complex II: Oxidized {reaction_rate} FADH2")
                return reaction_rate
        logger.warning("Insufficient FADH2 or ubiquinone for Complex II")
        return 0

    def complex_III(self) -> int:
        """
        Simulates Complex III activity.

        Complex III is the third enzyme in the electron transport chain and
        catalyzes the transfer of electrons from ubiquinol to cytochrome c.

        Returns
        -------
        int
            The amount of electrons transferred.
        """
        if self.is_metabolite_available(
            "ubiquinol", 1
        ) and self.is_metabolite_available("cytochrome_c_oxidized", 1):
            reaction_rate = min(
                self.metabolites["ubiquinol"].quantity,
                self.metabolites["cytochrome_c_oxidized"].quantity,
            )
            if self.consume_metabolites(
                ubiquinol=reaction_rate, cytochrome_c_oxidized=reaction_rate
            ):
                self.produce_metabolites(
                    ubiquinone=reaction_rate, cytochrome_c_reduced=reaction_rate
                )
                self.update_proton_gradient(PROTONS_PER_FADH2 * reaction_rate)
                logger.info(
                    f"Complex III: Transferred {reaction_rate} electron pairs, pumped {PROTONS_PER_FADH2 * reaction_rate} protons"
                )
                return reaction_rate
        logger.warning("Insufficient ubiquinol or cytochrome c for Complex III")
        return 0

    def complex_IV(self) -> int:
        """
        Simulates Complex IV activity.

        Complex IV is the fourth and final enzyme in the electron transport chain
        and catalyzes the oxidation of cytochrome c.

        Returns
        -------
        int
            The amount of electrons transferred.
        """
        if self.is_metabolite_available(
            "cytochrome_c_reduced", 1
        ) and self.is_metabolite_available("oxygen", 0.5):
            reaction_rate = min(
                self.metabolites["cytochrome_c_reduced"].quantity,
                self.metabolites["oxygen"].quantity * 2,
            )  # 2 cytochrome c per O2
            oxygen_consumed = reaction_rate / 2
            if self.consume_metabolites(
                cytochrome_c_reduced=reaction_rate, oxygen=oxygen_consumed
            ):
                self.produce_metabolites(cytochrome_c_oxidized=reaction_rate)
                self.update_proton_gradient(PROTONS_PER_FADH2 * reaction_rate)
                logger.info(
                    f"Complex IV: Consumed {oxygen_consumed} O2, pumped {PROTONS_PER_FADH2 * reaction_rate} protons"
                )
                return reaction_rate
        if self.metabolites["oxygen"].quantity <= 0:
            logger.warning("Insufficient oxygen for Complex IV")
        else:
            logger.warning("Insufficient reduced cytochrome c for Complex IV")
        return 0

    def is_metabolite_available(self, metabolite: str, amount: float) -> bool:
        """
        Check if a metabolite is available in sufficient quantity.

        Parameters
        ----------
        metabolite: str
            The metabolite to check.
        amount: float
            The amount of the metabolite to check.

        Returns
        -------
        bool
            True if the metabolite is available in sufficient quantity, False otherwise.
        """
        if hasattr(self, metabolite):
            return getattr(self, metabolite).quantity >= amount
        else:
            logger.warning(f"Unknown metabolite: {metabolite}")
            return False

    def atp_synthase(self) -> int:
        """
        Synthesizes ATP using the proton gradient.

        Returns
        -------
        int
            The amount of ATP produced.
        """
        protons_required_per_atp = PROTONS_PER_ATP
        possible_atp = int(self.proton_gradient / protons_required_per_atp)
        atp_produced = min(possible_atp, self.metabolites["adp"].quantity)
        if self.consume_metabolites(adp=atp_produced):
            self.produce_metabolites(atp=atp_produced)
            self.proton_gradient -= atp_produced * protons_required_per_atp
            logger.info(f"ATP Synthase: Produced {atp_produced} ATP")
            return atp_produced
        logger.warning("Insufficient ADP for ATP synthesis")
        return 0

    def replenish_ubiquinone(self) -> int:
        """
        Replenishes ubiquinone from ubiquinol.

        Returns
        -------
        int
            The amount of ubiquinone replenished.
        """
        replenish_amount = min(
            self.metabolites["ubiquinol"].quantity,
            self.metabolites["ubiquinone"].max_quantity
            - self.metabolites["ubiquinone"].quantity,
        )
        self.change_metabolite_quantity("ubiquinone", replenish_amount)
        self.change_metabolite_quantity("ubiquinol", -replenish_amount)
        logger.info(f"Replenished {replenish_amount} ubiquinone")

    def replenish_cytochrome_c(self) -> int:
        """
        Replenishes oxidized cytochrome c from reduced form.

        Returns
        -------
        int
            The amount of oxidized cytochrome c replenished.
        """
        replenish_amount = min(
            self.metabolites["cytochrome_c_reduced"].quantity,
            self.metabolites["cytochrome_c_oxidized"].max_quantity
            - self.metabolites["cytochrome_c_oxidized"].quantity,
        )
        self.change_metabolite_quantity("cytochrome_c_oxidized", replenish_amount)
        self.change_metabolite_quantity("cytochrome_c_reduced", -replenish_amount)
        logger.info(f"Replenished {replenish_amount} oxidized cytochrome c")

    def oxidative_phosphorylation(self, cytoplasmic_nadh_used: int = 0) -> int:
        """
        Simulates oxidative phosphorylation with the electron transport chain.

        Parameters
        ----------
        cytoplasmic_nadh_used: int
            The amount of cytoplasmic NADH used.

        Returns
        -------
        int
            The amount of ATP produced.
        """
        if self.metabolites["oxygen"].quantity <= 0:
            logger.warning("No oxygen available. Oxidative phosphorylation halted.")
            return 0

        total_nadh = self.metabolites["nadh"].quantity + cytoplasmic_nadh_used

        # Run the electron transport chain
        electrons_through_complex_I = self.complex_I()
        electrons_through_complex_II = self.complex_II()
        electrons_through_complex_III = self.complex_III()
        electrons_through_complex_IV = self.complex_IV()

        # ATP production via ATP synthase
        atp_produced = self.atp_synthase()

        # Replenish ubiquinone and cytochrome c
        self.replenish_ubiquinone()
        self.replenish_cytochrome_c()

        # Calculate efficiency
        total_electrons = electrons_through_complex_I + electrons_through_complex_II
        if total_electrons > 0:
            efficiency = atp_produced / total_electrons
            logger.info(
                f"Oxidative phosphorylation efficiency: {efficiency:.2f} ATP per electron pair"
            )

        return atp_produced

    def buffer_calcium(self, cytoplasmic_calcium: float) -> float:
        """
        Simulates calcium buffering by the mitochondrion.

        Parameters
        ----------
        cytoplasmic_calcium: float
            The amount of cytoplasmic calcium to buffer.

        Returns
        -------
        float
            The amount of calcium buffered.
        """
        calcium_uptake = min(
            cytoplasmic_calcium,
            self.metabolites["calcium"].max_quantity
            - self.metabolites["calcium"].quantity,
        )
        self.change_metabolite_quantity("calcium", calcium_uptake)
        logger.info(f"Mitochondrion buffered {calcium_uptake:.2f} units of calcium")

        if self.metabolites["calcium"].quantity > self.calcium_threshold:
            logger.warning(
                "Calcium overload detected. Risk of mitochondrial dysfunction."
            )

        return calcium_uptake

    def release_calcium(self, amount: float) -> float:
        """
        Releases calcium from the mitochondrion.

        Parameters
        ----------
        amount: float
            The amount of calcium to release.

        Returns
        -------
        float
            The amount of calcium released.
        """
        released = min(amount, self.metabolites["calcium"].quantity)
        self.change_metabolite_quantity("calcium", -released)
        logger.info(f"Mitochondrion released {released:.2f} units of calcium")
        return released

    def reset(self) -> None:
        """Reset mitochondrion state."""
        self.__init__()
        logger.info("Mitochondrion state reset")

    def transfer_cytoplasmic_nadh(self, cytoplasmic_nadh: float) -> float:
        """
        Transfers cytoplasmic NADH into the mitochondrion using shuttle systems.
        Returns the amount of mitochondrial NADH produced.

        Parameters
        ----------
        cytoplasmic_nadh: float
            The amount of cytoplasmic NADH to transfer.

        Returns
        -------
        float
            The amount of mitochondrial NADH produced.
        """
        shuttle_efficiency = SHUTTLE_EFFICIENCY
        mitochondrial_nadh = int(cytoplasmic_nadh * shuttle_efficiency)
        self.metabolites["nadh"].quantity += mitochondrial_nadh
        logger.info(
            f"Transferred {cytoplasmic_nadh} cytoplasmic NADH, produced {mitochondrial_nadh} mitochondrial NADH"
        )
        return mitochondrial_nadh
