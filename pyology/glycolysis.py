import logging
import math
from typing import TYPE_CHECKING

from pyology.common_reactions import GlycolysisReactions
from pyology.reaction import Reaction

from .exceptions import GlycolysisError, MetaboliteError, ReactionError
from .pathway import Pathway

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .organelle import Organelle


class Glycolysis(Pathway):
    """
    Class representing the glycolysis pathway.
    """

    time_step = 0.1  # Default time step in seconds

    reactions = GlycolysisReactions

    @classmethod
    def perform(
        cls, organelle: "Organelle", glucose_units: float
    ) -> tuple[float, float]:
        """
        Perform glycolysis on the given number of glucose units.

        Glycolysis is the process by which glucose is converted into pyruvate,
        releasing energy in the form of ATP.

        Parameters
        ----------
        organelle : Organelle
            The organelle where glycolysis takes place.
        glucose_units : float
            The number of glucose units to process.

        Returns
        -------
        tuple[float, float]
            The amount of net ATP produced and pyruvate produced.

        Raises
        ------
        GlycolysisError
            If glycolysis fails at any step or if the stoichiometry is incorrect.
        """
        logger.info(f"Starting glycolysis with {glucose_units} glucose units")
        logger.info(f"Initial metabolite levels: {organelle.metabolites.quantities}")
        glucose_units = math.floor(glucose_units)
        if glucose_units <= 0:
            raise GlycolysisError("The number of glucose units must be positive.")

        try:
            initial_atp = organelle.get_metabolite_quantity("ATP")
            initial_adp = organelle.get_metabolite_quantity("ADP")
            initial_amp = organelle.get_metabolite_quantity("AMP")
            initial_total = initial_atp + initial_adp + initial_amp
            logger.info(
                f"Initial ATP: {initial_atp}, Initial ADP: {initial_adp}, Initial AMP: {initial_amp}"
            )

            # Investment phase
            cls.investment_phase(organelle, glucose_units)

            # Yield phase
            cls.yield_phase(organelle, glucose_units * 2)  # 2 G3P per glucose

            # Add NAD+ regeneration step
            cls.regenerate_nad(organelle)

            final_atp = organelle.get_metabolite_quantity("ATP")
            final_adp = organelle.get_metabolite_quantity("ADP")
            final_amp = organelle.get_metabolite_quantity("AMP")
            final_total = final_atp + final_adp + final_amp
            logger.info(
                f"Final ATP: {final_atp}, Final ADP: {final_adp}, Final AMP: {final_amp}"
            )

            net_atp_produced = final_atp - initial_atp
            logger.info(f"Net ATP produced: {net_atp_produced}")

            # Ensure conservation of adenine nucleotides
            if abs(final_total - initial_total) > 1e-6:
                logger.warning(
                    f"Adenine nucleotide imbalance detected. Initial: {initial_total}, Final: {final_total}"
                )
                # Adjust ATP and ADP to maintain balance
                excess = final_total - initial_total
                atp_adjustment = min(excess, final_atp - initial_atp)
                organelle.metabolites["ATP"].quantity -= atp_adjustment
                adp_adjustment = excess - atp_adjustment
                if adp_adjustment > 0:
                    organelle.metabolites["ADP"].quantity -= adp_adjustment
                else:
                    organelle.metabolites["ADP"].quantity += abs(adp_adjustment)
                logger.info(
                    f"Adjusted ATP by -{atp_adjustment} and ADP by {-adp_adjustment} to maintain adenine nucleotide balance"
                )

            # Recalculate net ATP produced after adjustments
            final_atp = organelle.get_metabolite_quantity("ATP")
            net_atp_produced = final_atp - initial_atp
            logger.info(f"Adjusted net ATP produced: {net_atp_produced}")

            # Calculate pyruvate produced (2 pyruvate per glucose)
            pyruvate_produced = 2 * glucose_units

            logger.info(f"Glycolysis completed. Produced {pyruvate_produced} pyruvate.")
            logger.info(f"Final metabolite levels: {organelle.metabolites.quantities}")

            return net_atp_produced, pyruvate_produced

        except Exception as e:
            logger.error(f"Error during glycolysis: {str(e)}")
            raise GlycolysisError(f"Glycolysis failed: {str(e)}")

    @classmethod
    def investment_phase(cls, organelle, glucose_units):
        """
        Perform the investment phase of glycolysis (steps 1-5) for multiple glucose units.
        """
        logger.info(f"Starting investment phase with {glucose_units} glucose units")
        initial_atp = organelle.get_metabolite_quantity("ATP")

        for i in range(glucose_units):
            logger.info(
                f"🔄🔄🔄 Processing glucose unit {i+1} of {glucose_units} 🔄🔄🔄"
            )
            try:
                # Steps 1-4 occur once per glucose molecule
                cls.reactions.hexokinase.execute(organelle=organelle)
                cls.reactions.phosphoglucose_isomerase.execute(organelle=organelle)
                cls.reactions.phosphofructokinase.execute(organelle=organelle)
                cls.reactions.aldolase.execute(organelle=organelle)

                # Step 5 occurs once to convert DHAP to G3P
                cls.reactions.triose_phosphate_isomerase.execute(organelle=organelle)

                current_atp = organelle.get_metabolite_quantity("ATP")
                logger.info(f"ATP after processing glucose unit {i+1}: {current_atp}")

            except ReactionError as e:
                logger.error(f"Investment phase failed at glucose unit {i+1}: {str(e)}")
                raise GlycolysisError(f"Investment phase failed: {str(e)}")

        final_atp = organelle.get_metabolite_quantity("ATP")
        logger.info(f"ATP consumed in investment phase: {initial_atp - final_atp}")

    @classmethod
    def yield_phase(cls, organelle, g3p_units):
        """
        Perform the yield phase of glycolysis (steps 6-10) for multiple G3P units.
        """
        logger.info(f"Starting yield phase with {g3p_units} G3P units")
        initial_atp = organelle.get_metabolite_quantity("ATP")

        for i in range(g3p_units):
            logger.info(f"🍀🍀🍀 Processing G3P unit {i+1} of {g3p_units} 🍀🍀🍀")
            try:
                cls.reactions.glyceraldehyde_3_phosphate_dehydrogenase.execute(
                    organelle=organelle
                )
                cls.reactions.phosphoglycerate_kinase.execute(organelle=organelle)
                cls.reactions.phosphoglycerate_mutate.execute(
                    organelle=organelle
                )  # Changed from phosphoglycerate_mutate
                cls.reactions.enolase.execute(organelle=organelle)
                cls.reactions.pyruvate_kinase.execute(organelle=organelle)

                current_atp = organelle.get_metabolite_quantity("ATP")
                logger.info(f"ATP after processing G3P unit {i+1}: {current_atp}")

            except ReactionError as e:
                raise GlycolysisError(f"Yield phase failed at G3P unit {i+1}: {str(e)}")

        final_atp = organelle.get_metabolite_quantity("ATP")
        logger.info(f"ATP produced in yield phase: {final_atp - initial_atp}")

    @classmethod
    def enolase_reaction(cls, organelle: "Organelle"):
        initial_atp = organelle.get_metabolite_quantity("ATP")
        enolase = GlycolysisReactions.enolase
        phosphoglycerate_2 = organelle.get_metabolite_quantity("phosphoglycerate_2")

        while phosphoglycerate_2 > 0:
            reaction_result = enolase.execute(organelle)
            if reaction_result == 0:
                break
            phosphoglycerate_2 = organelle.get_metabolite_quantity("phosphoglycerate_2")

        final_atp = organelle.get_metabolite_quantity("ATP")
        logger.info(f"ATP change in enolase reaction: {final_atp - initial_atp}")
        logger.info(
            f"Enolase reaction completed. Remaining 2-phosphoglycerate: {phosphoglycerate_2}"
        )

    @classmethod
    def phosphoglycerate_kinase(cls, organelle):
        """1,3-Bisphosphoglycerate to 3-Phosphoglycerate"""
        reaction = cls.reactions.phosphoglycerate_kinase
        reaction.execute(organelle)

    @classmethod
    def pyruvate_kinase(cls, organelle):
        """Phosphoenolpyruvate to Pyruvate"""
        reaction = cls.reactions.pyruvate_kinase
        reaction.execute(organelle)

    @classmethod
    def regenerate_nad(cls, organelle):
        """Regenerate NAD+ via Lactate Dehydrogenase reaction"""
        nadh_quantity = organelle.get_metabolite_quantity("NADH")
        pyruvate_quantity = organelle.get_metabolite_quantity("pyruvate")

        reaction_amount = min(nadh_quantity, pyruvate_quantity)

        if reaction_amount > 0:
            cls.reactions.lactate_dehydrogenase.execute(organelle)
            logger.info(f"Regenerated {reaction_amount} NAD+ via Lactate Dehydrogenase")


def create_glycolysis_reactions():
    reactions = [
        # ... (other reactions remain unchanged)
        Reaction(
            name="Phosphoglycerate kinase",
            substrates={"1,3-bisphosphoglycerate": 1, "ADP": 1},
            products={"3-phosphoglycerate": 1, "ATP": 1},
            rate=1.0,
        ),
        # ... (other reactions remain unchanged)
        Reaction(
            name="Pyruvate kinase",
            substrates={"phosphoenolpyruvate": 1, "ADP": 1},
            products={"pyruvate": 1, "ATP": 1},
            rate=1.0,
        ),
    ]
    return reactions
