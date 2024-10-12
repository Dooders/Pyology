import logging
import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List

from constants import *

logger = logging.getLogger(__name__)


class Organelle:
    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.metabolites: Dict[str, Metabolite] = {}

    def add_metabolite(self, name: str, quantity: int, max_quantity: int):
        self.metabolites[name] = Metabolite(name, quantity, max_quantity)

    def change_metabolite_quantity(self, metabolite_name: str, amount: float):
        if metabolite_name in self.metabolites:
            metabolite = self.metabolites[metabolite_name]
            metabolite.quantity = max(
                0, min(metabolite.quantity + amount, metabolite.max_quantity)
            )
        else:
            self.logger.warning(f"Unknown metabolite: {metabolite_name}")

    def is_metabolite_available(self, metabolite: str, amount: float) -> bool:
        if metabolite in self.metabolites:
            return self.metabolites[metabolite].quantity >= amount
        else:
            self.logger.warning(f"Unknown metabolite: {metabolite}")
            return False

    def consume_metabolites(self, **metabolites: Dict[str, float]):
        for metabolite, amount in metabolites.items():
            if self.is_metabolite_available(metabolite, amount):
                self.change_metabolite_quantity(metabolite, -amount)
            else:
                self.logger.warning(f"Insufficient {metabolite} for reaction")
                return False
        return True

    def produce_metabolites(self, **metabolites: Dict[str, float]):
        for metabolite, amount in metabolites.items():
            self.change_metabolite_quantity(metabolite, amount)


@dataclass
class Enzyme:
    name: str
    activity: float = field(default=1.0)
    vmax: float = field(default=1.0)
    km: float = field(default=0.1)


@dataclass
class Metabolite:
    name: str = field(default="")
    quantity: int = field(default=0)
    max_quantity: int = field(default=1000)


@dataclass
class Effector:
    name: str = field(default="")
    concentration: float = field(default=0)
    Ki: float = field(default=0.1)  # Inhibition constant
    Ka: float = field(default=0.1)  # Activation constant


def michaelis_menten(substrate_conc: float, vmax: float, km: float) -> float:
    """Calculates reaction rate using the Michaelis-Menten equation."""
    return vmax * substrate_conc / (km + substrate_conc)


def allosteric_regulation(
    base_activity: float, inhibitors: List[Effector], activators: List[Effector]
) -> float:
    """Calculates enzyme activity considering inhibitors and activators."""
    inhibition_factor = 1
    for inhibitor in inhibitors:
        inhibition_factor *= 1 / (1 + inhibitor.concentration / inhibitor.Ki)
    activation_factor = 1
    for activator in activators:
        activation_factor *= 1 + activator.concentration / activator.Ka
    return base_activity * inhibition_factor * activation_factor


def hill_equation(substrate_conc: float, Vmax: float, K: float, n: float) -> float:
    """Calculates reaction rate using the Hill equation for cooperative binding."""
    return Vmax * (substrate_conc**n) / (K**n + substrate_conc**n)


class GlycolysisSteps(Enum):
    STEP1_HEXOKINASE = "Hexokinase"
    STEP2_PHOSPHOGLUCOSE_ISOMERASE = "Phosphoglucose Isomerase"
    STEP3_PHOSPHOFRUCTOKINASE = "Phosphofructokinase"
    STEP4_ALDOLASE = "Aldolase"
    STEP5_TRIOSE_PHOSPHATE_ISOMERASE = "Triose Phosphate Isomerase"
    STEP6_GLYCERALDEHYDE_3_PHOSPHATE_DEHYDROGENASE = (
        "Glyceraldehyde 3-Phosphate Dehydrogenase"
    )
    STEP7_PHOSPHOGLYCERATE_KINASE = "Phosphoglycerate Kinase"
    STEP8_PHOSPHOGLYCERATE_MUTASE = "Phosphoglycerate Mutase"
    STEP9_ENOLASE = "Enolase"
    STEP10_PYRUVATE_KINASE = "Pyruvate Kinase"


class Cytoplasm(Organelle):
    def __init__(self):
        super().__init__()
        self.add_metabolite("glucose", 0, 1000)
        self.add_metabolite("atp", 0, 1000)
        self.add_metabolite("adp", 0, 1000)
        self.add_metabolite("nad", 10, 1000)  # Starting NAD+ molecules
        self.add_metabolite("nadh", 0, 1000)
        self.add_metabolite("pyruvate", 0, 1000)
        self.glycolysis_rate = 1.0  # Base glycolysis rate

    def glycolysis(self, glucose_units):
        """
        Models glycolysis in a step-wise manner, including ATP and NADH
        production.
        """
        self.change_metabolite_quantity("glucose", glucose_units)
        logger.info(
            f"Starting glycolysis with {self.get_metabolite_quantity('glucose')} units of glucose"
        )

        for step in GlycolysisSteps:
            getattr(self, f"{step.name.lower()}")()

        logger.info(
            f"Glycolysis complete. Produced {self.get_metabolite_quantity('pyruvate')} pyruvate molecules"
        )
        return self.get_metabolite_quantity("pyruvate")

    def is_metabolite_available(self, metabolite: str, amount: float) -> bool:
        """Check if a metabolite is available in sufficient quantity."""
        if metabolite in self.metabolites:
            return self.metabolites[metabolite].quantity >= amount
        else:
            logger.warning(f"Unknown metabolite: {metabolite}")
            return False

    def consume_metabolites(self, **metabolites: Dict[str, float]):
        """Consume multiple metabolites at once."""
        for metabolite, amount in metabolites.items():
            if self.is_metabolite_available(metabolite, amount):
                self.metabolites[metabolite].quantity -= amount
            else:
                logger.warning(f"Insufficient {metabolite} for reaction")
                return False
        return True

    def produce_metabolites(self, **metabolites: Dict[str, float]):
        """Produce multiple metabolites at once."""
        for metabolite, amount in metabolites.items():
            self.metabolites[metabolite].quantity += amount

    def step1_hexokinase(self):
        if self.consume_metabolites(glucose=1, atp=1):
            self.produce_metabolites(adp=1)
            logger.info(
                f"Step 1: {GlycolysisSteps.STEP1_HEXOKINASE.value} - Glucose phosphorylation"
            )

    def step2_phosphoglucose_isomerase(self):
        logger.info(
            f"Step 2: {GlycolysisSteps.STEP2_PHOSPHOGLUCOSE_ISOMERASE.value} - Isomerization"
        )

    def step3_phosphofructokinase(self):
        if self.consume_metabolites(atp=1):
            self.produce_metabolites(adp=1)
            logger.info(
                f"Step 3: {GlycolysisSteps.STEP3_PHOSPHOFRUCTOKINASE.value} - Phosphorylation"
            )

    def step4_aldolase(self):
        logger.info(
            f"Step 4: {GlycolysisSteps.STEP4_ALDOLASE.value} - Splitting fructose-1,6-bisphosphate"
        )

    def step5_triose_phosphate_isomerase(self):
        logger.info(
            f"Step 5: {GlycolysisSteps.STEP5_TRIOSE_PHOSPHATE_ISOMERASE.value} - Isomerization"
        )

    def step6_glyceraldehyde_3_phosphate_dehydrogenase(self):
        if self.consume_metabolites(nad=2):
            self.produce_metabolites(nadh=2)
            logger.info(
                f"Step 6: {GlycolysisSteps.STEP6_GLYCERALDEHYDE_3_PHOSPHATE_DEHYDROGENASE.value} - Oxidation and phosphorylation"
            )

    def step7_phosphoglycerate_kinase(self):
        if self.consume_metabolites(adp=2):
            self.produce_metabolites(atp=2)
            logger.info(
                f"Step 7: {GlycolysisSteps.STEP7_PHOSPHOGLYCERATE_KINASE.value} - ATP generation"
            )

    def step8_phosphoglycerate_mutase(self):
        logger.info(
            f"Step 8: {GlycolysisSteps.STEP8_PHOSPHOGLYCERATE_MUTASE.value} - Shifting phosphate group"
        )

    def step9_enolase(self):
        logger.info(f"Step 9: {GlycolysisSteps.STEP9_ENOLASE.value} - Dehydration")

    def step10_pyruvate_kinase(self):
        if self.consume_metabolites(adp=2):
            self.produce_metabolites(atp=2, pyruvate=2)
            logger.info(
                f"Step 10: {GlycolysisSteps.STEP10_PYRUVATE_KINASE.value} - ATP generation and pyruvate formation"
            )

    def reset(self):
        self.__init__()
        logger.info("Cytoplasm state reset")

    def get_metabolite_quantity(self, metabolite_name: str) -> int:
        """Get the quantity of a specific metabolite."""
        if metabolite_name in self.metabolites:
            return self.metabolites[metabolite_name].quantity
        else:
            logger.warning(f"Unknown metabolite: {metabolite_name}")
            return 0


class Mitochondrion(Organelle):
    """
    Simulates the electron transport chain with individual complexes.

    Proton gradient and ATP synthesis are modeled well, including a proton leak
    for added realism.

    Calcium buffering and its effect on mitochondrial function are included,
    which adds another layer of detail.

    Feedback inhibition for ATP levels is implemented, mimicking real cellular
    regulation mechanisms.
    """

    def __init__(self):
        super().__init__()
        self.add_metabolite("nadh", 0, 1000)
        self.add_metabolite("fadh2", 0, 1000)
        self.add_metabolite("atp", 0, 1000)
        self.add_metabolite("adp", 100, 1000)
        self.add_metabolite("oxygen", 1000, 1000)
        self.add_metabolite("ubiquinone", 100, 1000)
        self.add_metabolite("ubiquinol", 0, 1000)
        self.add_metabolite("cytochrome_c_oxidized", 100, 1000)
        self.add_metabolite("cytochrome_c_reduced", 0, 1000)
        self.add_metabolite("co2", 0, 1000000)
        self.add_metabolite("calcium", 0, 1000)
        self.add_metabolite("oxaloacetate", 0, 1000)

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

    def change_metabolite_quantity(self, metabolite_name: str, amount: float):
        """Centralized method to change metabolite quantities."""
        if metabolite_name in self.metabolites:
            metabolite = self.metabolites[metabolite_name]
            metabolite.quantity = max(
                0, min(metabolite.quantity + amount, metabolite.max_quantity)
            )
        else:
            logger.warning(f"Unknown metabolite: {metabolite_name}")

    def consume_metabolites(self, **metabolites: Dict[str, float]):
        """Consume multiple metabolites at once."""
        for metabolite, amount in metabolites.items():
            if self.is_metabolite_available(metabolite, amount):
                self.metabolites[metabolite].quantity -= amount
            else:
                logger.warning(f"Insufficient {metabolite} for reaction")
                return False
        return True

    def produce_metabolites(self, **metabolites: Dict[str, float]):
        """Produce multiple metabolites at once."""
        for metabolite, amount in metabolites.items():
            self.metabolites[metabolite].quantity += amount

    def krebs_cycle_process(self, acetyl_coa_amount: int):
        """Processes acetyl-CoA through the Krebs cycle"""
        logger.info(
            f"Processing {acetyl_coa_amount} units of acetyl-CoA through the Krebs cycle"
        )

        self.krebs_cycle.add_substrate("acetyl_coa", acetyl_coa_amount)

        # Ensure there's enough oxaloacetate to start the cycle
        if self.metabolites["oxaloacetate"].quantity < acetyl_coa_amount:
            self.krebs_cycle.add_substrate(
                "oxaloacetate",
                acetyl_coa_amount - self.metabolites["oxaloacetate"].quantity,
            )

        for _ in range(acetyl_coa_amount):
            self.krebs_cycle.run_cycle()

        # Transfer the products to the mitochondrion
        self.change_metabolite_quantity("nadh", self.krebs_cycle.cofactors["nadh"])
        self.change_metabolite_quantity("fadh2", self.krebs_cycle.cofactors["fadh2"])
        self.change_metabolite_quantity("atp", self.krebs_cycle.cofactors["gtp"])

        # Reset the Krebs cycle for the next round
        self.krebs_cycle.reset()

        return self.krebs_cycle.cofactors["nadh"] + self.krebs_cycle.cofactors["fadh2"]

    def pyruvate_to_acetyl_coa(self, pyruvate_amount: int) -> int:
        """Converts pyruvate to acetyl-CoA."""
        logger.info(f"Converting {pyruvate_amount} units of pyruvate to acetyl-CoA")
        acetyl_coa_produced = pyruvate_amount
        self.change_metabolite_quantity("nadh", pyruvate_amount)
        self.change_metabolite_quantity("co2", pyruvate_amount)
        return acetyl_coa_produced

    def cellular_respiration(self, pyruvate_amount: int):
        """Simulates the entire cellular respiration process with feedback inhibition"""
        if self.metabolites["oxygen"].quantity <= 0:
            logger.warning("No oxygen available. Cellular respiration halted.")
            return 0

        acetyl_coa = self.pyruvate_to_acetyl_coa(pyruvate_amount)

        # Implement feedback inhibition
        atp_inhibition_factor = 1 / (
            1 + self.metabolites["atp"].quantity / 1000
        )  # Example threshold
        self.krebs_cycle.enzymes["citrate_synthase"].activity *= atp_inhibition_factor
        self.krebs_cycle.enzymes[
            "isocitrate_dehydrogenase"
        ].activity *= atp_inhibition_factor

        krebs_products = self.krebs_cycle_process(acetyl_coa)

        # Transfer NADH and FADH2 from Krebs cycle to ETC
        self.change_metabolite_quantity("nadh", self.krebs_cycle.cofactors["nadh"])
        self.change_metabolite_quantity("fadh2", self.krebs_cycle.cofactors["fadh2"])

        # Check ADP availability
        if self.metabolites["adp"].quantity < 10:  # Arbitrary threshold
            logger.warning("Low ADP levels. Oxidative phosphorylation may be limited.")

        atp_produced = self.oxidative_phosphorylation()

        # Add ATP from substrate-level phosphorylation in Krebs cycle
        atp_produced += self.krebs_cycle.cofactors["gtp"]  # GTP is equivalent to ATP

        return atp_produced

    def calculate_proton_leak(self):
        """Calculate the proton leak using a logistic function."""
        relative_gradient = self.proton_gradient / self.max_proton_gradient
        leak = self.leak_rate / (
            1
            + math.exp(
                -self.leak_steepness * (self.proton_gradient - self.leak_midpoint)
            )
        )
        return leak

    def update_proton_gradient(self, protons_pumped):
        """Update the proton gradient considering nonlinear leak."""
        self.proton_gradient += protons_pumped
        leak = self.calculate_proton_leak()
        self.proton_gradient = max(0, self.proton_gradient - leak)
        logger.info(f"Proton gradient: {self.proton_gradient:.2f}, Leak: {leak:.2f}")

    def complex_I(self):
        """Simulates Complex I activity."""
        if self.is_metabolite_available("nadh", 1) and self.is_metabolite_available(
            "ubiquinone", 1
        ):
            reaction_rate = min(
                self.metabolites["nadh"].quantity,
                self.metabolites["ubiquinone"].quantity,
            )
            if self.consume_metabolites(nadh=reaction_rate, ubiquinone=reaction_rate):
                self.produce_metabolites(ubiquinol=reaction_rate)
                self.proton_gradient += PROTONS_PER_NADH * reaction_rate
                logger.info(
                    f"Complex I: Oxidized {reaction_rate} NADH, pumped {PROTONS_PER_NADH * reaction_rate} protons"
                )
                return reaction_rate
        logger.warning("Insufficient NADH or ubiquinone for Complex I")
        return 0

    def complex_II(self):
        """Simulates Complex II activity."""
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

    def complex_III(self):
        """Simulates Complex III activity."""
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
                self.proton_gradient += PROTONS_PER_FADH2 * reaction_rate
                logger.info(
                    f"Complex III: Transferred {reaction_rate} electron pairs, pumped {PROTONS_PER_FADH2 * reaction_rate} protons"
                )
                return reaction_rate
        logger.warning("Insufficient ubiquinol or cytochrome c for Complex III")
        return 0

    def complex_IV(self):
        """Simulates Complex IV activity."""
        if self.is_metabolite_available(
            "cytochrome_c_reduced", 1
        ) and self.is_metabolite_available("oxygen", 0.5):
            reaction_rate = min(
                self.metabolites["cytochrome_c_reduced"].quantity,
                self.metabolites["oxygen"].quantity * 2,
            )  # 2 cytochrome c per O2
            oxygen_consumed = reaction_rate // 2
            if self.consume_metabolites(
                cytochrome_c_reduced=reaction_rate, oxygen=oxygen_consumed
            ):
                self.produce_metabolites(cytochrome_c_oxidized=reaction_rate)
                self.proton_gradient += PROTONS_PER_FADH2 * reaction_rate
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
        """Check if a metabolite is available in sufficient quantity."""
        if hasattr(self, metabolite):
            return getattr(self, metabolite).quantity >= amount
        else:
            logger.warning(f"Unknown metabolite: {metabolite}")
            return False

    def atp_synthase(self):
        """Synthesizes ATP using the proton gradient."""
        protons_required_per_atp = PROTONS_PER_ATP
        possible_atp = int(self.proton_gradient / protons_required_per_atp)
        atp_produced = min(possible_atp, self.metabolites["adp"].quantity)
        self.change_metabolite_quantity("atp", atp_produced)
        self.change_metabolite_quantity("adp", -atp_produced)
        self.proton_gradient -= atp_produced * protons_required_per_atp
        logger.info(f"ATP Synthase: Produced {atp_produced} ATP")
        return atp_produced

    def replenish_ubiquinone(self):
        """Replenishes ubiquinone from ubiquinol"""
        replenish_amount = min(
            self.metabolites["ubiquinol"].quantity,
            self.metabolites["ubiquinone"].max_quantity
            - self.metabolites["ubiquinone"].quantity,
        )
        self.change_metabolite_quantity("ubiquinone", replenish_amount)
        self.change_metabolite_quantity("ubiquinol", -replenish_amount)
        logger.info(f"Replenished {replenish_amount} ubiquinone")

    def replenish_cytochrome_c(self):
        """Replenishes oxidized cytochrome c from reduced form"""
        replenish_amount = min(
            self.metabolites["cytochrome_c_reduced"].quantity,
            self.metabolites["cytochrome_c_oxidized"].max_quantity
            - self.metabolites["cytochrome_c_oxidized"].quantity,
        )
        self.change_metabolite_quantity("cytochrome_c_oxidized", replenish_amount)
        self.change_metabolite_quantity("cytochrome_c_reduced", -replenish_amount)
        logger.info(f"Replenished {replenish_amount} oxidized cytochrome c")

    def oxidative_phosphorylation(self, cytoplasmic_nadh_used: int = 0):
        """Simulates oxidative phosphorylation with the electron transport chain."""
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

    def buffer_calcium(self, cytoplasmic_calcium: int):
        """Simulates calcium buffering by the mitochondrion."""
        calcium_uptake = min(
            cytoplasmic_calcium,
            self.metabolites["calcium"].max_quantity
            - self.metabolites["calcium"].quantity,
        )
        self.change_metabolite_quantity("calcium", calcium_uptake)
        logger.info(f"Mitochondrion buffered {calcium_uptake} units of calcium")

        if self.metabolites["calcium"].quantity > self.calcium_threshold:
            logger.warning(
                "Calcium overload detected. Risk of mitochondrial dysfunction."
            )

        return calcium_uptake

    def release_calcium(self, amount: int):
        """Releases calcium from the mitochondrion."""
        released = min(amount, self.metabolites["calcium"].quantity)
        self.change_metabolite_quantity("calcium", -released)
        logger.info(f"Mitochondrion released {released} units of calcium")
        return released

    def reset(self):
        """Reset mitochondrion state."""
        self.__init__()
        logger.info("Mitochondrion state reset")

    def transfer_cytoplasmic_nadh(self, cytoplasmic_nadh: int) -> int:
        """
        Transfers cytoplasmic NADH into the mitochondrion using shuttle systems.
        Returns the amount of mitochondrial NADH produced.
        """
        shuttle_efficiency = SHUTTLE_EFFICIENCY
        mitochondrial_nadh = int(cytoplasmic_nadh * shuttle_efficiency)
        self.change_metabolite_quantity("nadh", mitochondrial_nadh)
        logger.info(
            f"Transferred {cytoplasmic_nadh} cytoplasmic NADH, produced {mitochondrial_nadh} mitochondrial NADH"
        )
        return mitochondrial_nadh


class KrebsCycle(Organelle):
    """
    Krebs cycle is modeled step-by-step, with each enzyme's activity influenced
    by effectors and inhibitors.

    Uses Michaelis-Menten kinetics and Hill equations, which are common for
    modeling enzyme-catalyzed reactions.

    Includes regulation of enzyme activity by cofactors like ATP and ADP, which
    is a realistic representation of metabolic regulation.
    """

    def __init__(self):
        super().__init__()
        self.add_metabolite("Acetyl-CoA", 0, 1000)
        self.add_metabolite("Oxaloacetate", 0, 1000)
        self.add_metabolite("Citrate", 0, 1000)
        self.add_metabolite("Isocitrate", 0, 1000)
        self.add_metabolite("α-Ketoglutarate", 0, 1000)
        self.add_metabolite("Succinyl-CoA", 0, 1000)
        self.add_metabolite("Succinate", 0, 1000)
        self.add_metabolite("Fumarate", 0, 1000)
        self.add_metabolite("Malate", 0, 1000)

        self.cofactors = {
            "nad": INITIAL_NAD,
            "nadh": 0,
            "fad": 100,
            "fadh2": 0,
            "coenzyme_a": 100,
            "atp": 100,
            "adp": 0,
            "gtp": 0,
            "gdp": 0,
            "co2": 0,
        }
        self.enzymes = {
            "citrate_synthase": Enzyme("Citrate Synthase"),
            "aconitase": Enzyme("Aconitase"),
            "isocitrate_dehydrogenase": Enzyme("Isocitrate Dehydrogenase"),
            "alpha_ketoglutarate_dehydrogenase": Enzyme(
                "α-Ketoglutarate Dehydrogenase"
            ),
            "succinyl_coa_synthetase": Enzyme("Succinyl-CoA Synthetase"),
            "succinate_dehydrogenase": Enzyme("Succinate Dehydrogenase"),
            "fumarase": Enzyme("Fumarase"),
            "malate_dehydrogenase": Enzyme("Malate Dehydrogenase"),
        }

    def is_metabolite_available(self, metabolite: str, amount: float) -> bool:
        """Check if a metabolite is available in sufficient quantity."""
        if metabolite in self.metabolites:
            return self.metabolites[metabolite] >= amount
        elif metabolite in self.cofactors:
            return self.cofactors[metabolite] >= amount
        else:
            logger.warning(f"Unknown metabolite: {metabolite}")
            return False

    def consume_metabolites(self, **metabolites: Dict[str, float]):
        """Consume multiple metabolites at once."""
        for metabolite, amount in metabolites.items():
            if self.is_metabolite_available(metabolite, amount):
                if metabolite in self.metabolites:
                    self.metabolites[metabolite] -= amount
                elif metabolite in self.cofactors:
                    self.cofactors[metabolite] -= amount
            else:
                logger.warning(f"Insufficient {metabolite} for reaction")
                return False
        return True

    def produce_metabolites(self, **metabolites: Dict[str, float]):
        """Produce multiple metabolites at once."""
        for metabolite, amount in metabolites.items():
            if metabolite in self.metabolites:
                self.metabolites[metabolite] += amount
            elif metabolite in self.cofactors:
                self.cofactors[metabolite] += amount
            else:
                logger.warning(f"Unknown metabolite: {metabolite}")

    def step1_citrate_synthase(self):
        """Acetyl-CoA + Oxaloacetate to Citrate"""
        enzyme = self.enzymes["citrate_synthase"]
        substrate_conc = min(
            self.metabolites["acetyl_coa"], self.metabolites["oxaloacetate"]
        )
        reaction_rate = michaelis_menten(
            substrate_conc, enzyme.vmax * enzyme.activity, enzyme.km
        )

        if self.consume_metabolites(
            acetyl_coa=reaction_rate, oxaloacetate=reaction_rate
        ):
            self.produce_metabolites(citrate=reaction_rate, coenzyme_a=reaction_rate)
        else:
            logger.warning("Insufficient substrates for step 1")

    def step2_aconitase(self):
        """Citrate to Isocitrate"""
        enzyme = self.enzymes["aconitase"]
        substrate_conc = self.metabolites["citrate"]
        reaction_rate = michaelis_menten(
            substrate_conc, enzyme.vmax * enzyme.activity, enzyme.km
        )

        if self.consume_metabolites(citrate=reaction_rate):
            self.produce_metabolites(isocitrate=reaction_rate)
        else:
            logger.warning("Insufficient citrate for step 2")

    def step3_isocitrate_dehydrogenase(self):
        """Isocitrate to α-Ketoglutarate with allosteric regulation"""
        enzyme = self.enzymes["isocitrate_dehydrogenase"]
        substrate_conc = self.metabolites["isocitrate"]

        # Define effectors
        atp_effector = Effector("ATP", self.cofactors["atp"], Ki=100, Ka=1000)
        adp_effector = Effector("ADP", self.cofactors["adp"], Ki=1000, Ka=100)

        # Calculate regulated enzyme activity
        regulated_activity = allosteric_regulation(
            enzyme.activity,
            inhibitors=[atp_effector],
            activators=[adp_effector],
        )

        # Use Hill equation for cooperative binding
        n = 2  # Hill coefficient
        reaction_rate = hill_equation(
            substrate_conc, enzyme.vmax * regulated_activity, enzyme.km, n
        )

        if self.consume_metabolites(isocitrate=reaction_rate, nad=reaction_rate):
            self.produce_metabolites(
                α_ketoglutarate=reaction_rate, nadh=reaction_rate, co2=reaction_rate
            )
        else:
            logger.warning("Insufficient substrates or NAD⁺ for step 3")

    def step4_alpha_ketoglutarate_dehydrogenase(self):
        """α-Ketoglutarate to Succinyl-CoA"""
        enzyme = self.enzymes["alpha_ketoglutarate_dehydrogenase"]
        substrate_conc = self.metabolites["α-Ketoglutarate"]

        # Enzyme regulation
        atp_inhibition = self.cofactors["atp"] / 100
        nadh_inhibition = self.cofactors["nadh"] / 100
        succinyl_coa_inhibition = (
            self.metabolites["succinyl_coa"] / 10
        )  # Assuming max succinyl-CoA is 10
        enzyme_activity = (
            1 - (atp_inhibition + nadh_inhibition + succinyl_coa_inhibition) / 3
        )

        reaction_rate = michaelis_menten(
            substrate_conc,
            enzyme.vmax * enzyme_activity * enzyme.activity,
            enzyme.km,
        )

        if self.consume_metabolites(α_ketoglutarate=reaction_rate, nad=reaction_rate):
            self.produce_metabolites(
                succinyl_coa=reaction_rate, nadh=reaction_rate, co2=reaction_rate
            )
        else:
            logger.warning("Insufficient substrates or NAD⁺ for step 4")

    def step5_succinyl_coa_synthetase(self):
        """Succinyl-CoA to Succinate"""
        enzyme = self.enzymes["succinyl_coa_synthetase"]
        substrate_conc = self.metabolites["succinyl_coa"]
        reaction_rate = michaelis_menten(
            substrate_conc, enzyme.vmax * enzyme.activity, enzyme.km
        )

        if self.consume_metabolites(succinyl_coa=reaction_rate, gdp=reaction_rate):
            self.produce_metabolites(
                succinate=reaction_rate, gtp=reaction_rate, coenzyme_a=reaction_rate
            )
        else:
            logger.warning("Insufficient substrates or GDP for step 5")

    def step6_succinate_dehydrogenase(self):
        """Succinate to Fumarate"""
        enzyme = self.enzymes["succinate_dehydrogenase"]
        substrate_conc = self.metabolites["succinate"]
        reaction_rate = michaelis_menten(
            substrate_conc, enzyme.vmax * enzyme.activity, enzyme.km
        )

        if self.consume_metabolites(succinate=reaction_rate, fad=reaction_rate):
            self.produce_metabolites(fumarate=reaction_rate, fadh2=reaction_rate)
        else:
            logger.warning("Insufficient substrates or FAD for step 6")

    def step7_fumarase(self):
        """Fumarate to Malate"""
        enzyme = self.enzymes["fumarase"]
        substrate_conc = self.metabolites["fumarate"]
        reaction_rate = michaelis_menten(
            substrate_conc, enzyme.vmax * enzyme.activity, enzyme.km
        )

        if self.consume_metabolites(fumarate=reaction_rate):
            self.produce_metabolites(malate=reaction_rate)
        else:
            logger.warning("Insufficient fumarate for step 7")

    def step8_malate_dehydrogenase(self):
        """Malate to Oxaloacetate"""
        enzyme = self.enzymes["malate_dehydrogenase"]
        substrate_conc = self.metabolites["malate"]
        reaction_rate = michaelis_menten(
            substrate_conc, enzyme.vmax * enzyme.activity, enzyme.km
        )

        if self.consume_metabolites(malate=reaction_rate, nad=reaction_rate):
            self.produce_metabolites(oxaloacetate=reaction_rate, nadh=reaction_rate)
        else:
            logger.warning("Insufficient substrates or NAD⁺ for step 8")

    def run_cycle(self):
        self.step1_citrate_synthase()
        self.step2_aconitase()
        self.step3_isocitrate_dehydrogenase()
        self.step4_alpha_ketoglutarate_dehydrogenase()
        self.step5_succinyl_coa_synthetase()
        self.step6_succinate_dehydrogenase()
        self.step7_fumarase()
        self.step8_malate_dehydrogenase()

    def add_substrate(self, substrate: str, amount: float):
        """Add initial substrate to start the cycle"""
        if substrate in self.metabolites:
            self.metabolites[substrate] += amount
        elif substrate in self.cofactors:
            self.cofactors[substrate] += amount
        else:
            logger.warning(f"Unknown substrate: {substrate}")

    def display_state(self):
        """Display the current state of metabolites and cofactors"""
        print("Metabolites:")
        for metabolite, amount in self.metabolites.items():
            print(f"  {metabolite}: {amount:.2f}")
        print("\nCofactors:")
        for cofactor, amount in self.cofactors.items():
            print(f"  {cofactor}: {amount:.2f}")

    def reset(self):
        """Reset the Krebs cycle to its initial state"""
        self.__init__()


class Cell:
    def __init__(self):
        self.cytoplasm = Cytoplasm()
        self.mitochondrion = Mitochondrion()
        self.krebs_cycle = KrebsCycle()
        self.simulation_time = 0
        self.time_step = TIME_STEP
        self.cytoplasmic_calcium = Metabolite("Ca2+", 100, 1000)

    def produce_atp(self, glucose, simulation_duration=SIMULATION_DURATION):
        """Simulates ATP production in the entire cell over a specified duration."""
        initial_atp = (
            self.cytoplasm.metabolites["atp"].quantity
            + self.mitochondrion.metabolites["atp"].quantity
        )
        self.simulation_time = 0
        total_atp_produced = 0
        glucose_processed = 0

        while (
            glucose_processed < glucose and self.simulation_time < simulation_duration
        ):
            if self.mitochondrion.metabolites["oxygen"].quantity <= 0:
                logger.warning("Oxygen depleted. Stopping simulation.")
                break

            # Check ADP availability
            if (
                self.mitochondrion.metabolites["adp"].quantity < 10
            ):  # Arbitrary threshold
                logger.warning(
                    "Low ADP levels in mitochondrion. Transferring ADP from cytoplasm."
                )
                adp_transfer = min(
                    50, self.cytoplasm.metabolites["adp"].quantity
                )  # Transfer up to 50 ADP
                self.mitochondrion.metabolites["adp"].quantity += adp_transfer
                self.cytoplasm.metabolites["adp"].quantity -= adp_transfer

            # Implement feedback activation
            adp_activation_factor = (
                1 + self.cytoplasm.metabolites["adp"].quantity / 500
            )  # Example threshold
            self.cytoplasm.glycolysis_rate *= adp_activation_factor

            # Glycolysis with updated rate
            pyruvate = self.cytoplasm.glycolysis(1 * self.cytoplasm.glycolysis_rate)
            glucose_processed += 1 * self.cytoplasm.glycolysis_rate

            # Calculate ATP produced in glycolysis
            glycolysis_atp = (
                self.cytoplasm.metabolites["atp"].quantity
                - initial_atp
                + self.mitochondrion.metabolites["atp"].quantity
            )
            total_atp_produced += glycolysis_atp

            cytoplasmic_nadh = self.cytoplasm.metabolites["nadh"].quantity

            # NADH shuttle
            mitochondrial_nadh = self.mitochondrion.transfer_cytoplasmic_nadh(
                cytoplasmic_nadh
            )

            # Cellular respiration in mitochondrion
            mitochondrial_atp = self.mitochondrion.cellular_respiration(pyruvate)
            total_atp_produced += mitochondrial_atp

            # Transfer excess ATP from mitochondrion to cytoplasm
            atp_transfer = max(
                0, self.mitochondrion.metabolites["atp"].quantity - 100
            )  # Keep 100 ATP in mitochondrion
            self.cytoplasm.metabolites["atp"].quantity += atp_transfer
            self.mitochondrion.metabolites["atp"].quantity -= atp_transfer

            self.simulation_time += self.time_step

        atp_per_glucose = (
            total_atp_produced / glucose_processed if glucose_processed > 0 else 0
        )

        logger.info(
            f"Simulation completed. Time elapsed: {self.simulation_time:.2f} seconds"
        )
        logger.info(f"Glucose units processed: {glucose_processed}")
        logger.info(f"Total ATP produced: {total_atp_produced}")
        logger.info(f"ATP yield per glucose molecule: {atp_per_glucose:.2f}")
        logger.info(
            f"Remaining oxygen: {self.mitochondrion.metabolites['oxygen'].quantity}"
        )

        return total_atp_produced

    def reset(self):
        """Reset the entire cell state."""
        self.cytoplasm.reset()
        self.mitochondrion.reset()
        self.simulation_time = 0
        self.cytoplasmic_calcium = Metabolite(
            "Ca2+", 100, 1000
        )  # Reset cytoplasmic calcium
        logger.info("Cell state reset")


# Simulation code
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    cell = Cell()
    glucose_amounts = [1, 2, 5, 10]
    simulation_duration = SIMULATION_DURATION

    for glucose in glucose_amounts:
        logger.info(f"\nSimulating ATP production with {glucose} glucose units:")
        atp_produced = cell.produce_atp(glucose, simulation_duration)

        # Log the current state of the cell
        logger.info(f"Cytoplasm ATP: {cell.cytoplasm.metabolites['atp'].quantity}")
        logger.info(f"Cytoplasm NADH: {cell.cytoplasm.metabolites['nadh'].quantity}")
        logger.info(
            f"Mitochondrion ATP: {cell.mitochondrion.metabolites['atp'].quantity}"
        )
        logger.info(
            f"Mitochondrion NADH: {cell.mitochondrion.metabolites['nadh'].quantity}"
        )
        logger.info(
            f"Mitochondrion FADH2: {cell.mitochondrion.metabolites['fadh2'].quantity}"
        )
        logger.info(f"Simulation time: {cell.simulation_time:.2f} seconds")
        logger.info(
            f"Mitochondrial Calcium: {cell.mitochondrion.metabolites['calcium'].quantity}"
        )
        logger.info(f"Cytoplasmic Calcium: {cell.cytoplasmic_calcium.quantity}")
        logger.info(f"Proton Gradient: {cell.mitochondrion.proton_gradient:.2f}")

        # Reset the cell for the next simulation
        cell.reset()

    logger.info("Simulation complete.")
