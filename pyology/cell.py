import logging

from .constants import SIMULATION_DURATION, TIME_STEP
from .cytoplasm import Cytoplasm
from .mitochondrion import KrebsCycle, Mitochondrion
from .organelle import Organelle


class Cell(Organelle):
    name = "Cell"

    def __init__(self, logger=None, debug=False):
        super().__init__()
        self.atp = 300  # Initial ATP value
        self.cytoplasm = Cytoplasm(logger=logger, debug=debug)
        self.mitochondrion = Mitochondrion()
        self.krebs_cycle = KrebsCycle()
        self.simulation_time = 0
        self.time_step = TIME_STEP
        self.base_glycolysis_rate = 1.0  # This is now in the Cell class
        self.glycolysis_rate = self.base_glycolysis_rate  # Initialize glycolysis_rate
        self.initial_adenine_nucleotides = None
        self.debug = debug
        self.logger = logger

    def produce_atp(self, glucose, simulation_duration=SIMULATION_DURATION):
        #! Not being used, may not need
        """Simulates ATP production in the entire cell over a specified duration."""
        self.simulation_time = 0
        total_atp_produced = 0
        glucose_processed = 0

        while (
            glucose_processed < glucose and self.simulation_time < simulation_duration
        ):
            if self.metabolites["oxygen"].quantity <= 0:
                self.logger.warning("Oxygen depleted. Stopping simulation.")
                break

            # Calculate ATP at the start of the iteration
            atp_start = (
                self.metabolites["atp"].quantity + self.metabolites["atp"].quantity
            )

            # Check ADP availability
            if self.metabolites["adp"].quantity < 10:  # Arbitrary threshold
                self.logger.warning(
                    "Low ADP levels in mitochondrion. Transferring ADP from cytoplasm."
                )
                adp_transfer = min(
                    50, max(0, self.metabolites["adp"].quantity)
                )  # Transfer up to 50 ADP, but not less than 0
                self.metabolites["adp"].quantity += adp_transfer
                self.metabolites["adp"].quantity -= adp_transfer

            # Implement feedback activation
            adp_activation_factor = (
                1 + self.metabolites["adp"].quantity / 500
            )  # Example threshold
            self.glycolysis_rate = self.base_glycolysis_rate * adp_activation_factor

            # Glycolysis with updated rate
            glucose_consumed = min(
                1 * self.glycolysis_rate, self.metabolites["glucose"].quantity
            )
            pyruvate = self.glycolysis.perform_glycolysis(glucose_consumed)
            self.metabolites.add(
                "pyruvate", pyruvate, max_quantity=1000
            )  # Updated: added max_quantity
            self.logger.info(f"Transferred {pyruvate} pyruvate to mitochondrion")
            glucose_processed += glucose_consumed

            cytoplasmic_nadh = self.metabolites["nadh"].quantity

            # NADH shuttle
            mitochondrial_nadh = self.mitochondrion.transfer_cytoplasmic_nadh(
                cytoplasmic_nadh
            )

            # Cellular respiration in mitochondrion
            #! Pausing for now
            # mitochondrial_atp = self.mitochondrion.cellular_respiration(pyruvate)

            # Transfer excess ATP from mitochondrion to cytoplasm
            atp_transfer = max(
                0, self.metabolites["atp"].quantity - 100
            )  # Keep 100 ATP in mitochondrion
            self.metabolites["atp"].quantity += atp_transfer
            self.metabolites["atp"].quantity -= atp_transfer

            # Calculate ATP at the end of the iteration
            atp_end = (
                self.metabolites["atp"].quantity + self.metabolites["atp"].quantity
            )

            # Calculate ATP produced in this iteration
            delta_atp = atp_end - atp_start
            total_atp_produced += delta_atp

            self.simulation_time += self.time_step

        # ... existing code for logging results ...

        return total_atp_produced

    def produce_atp_generator(self, glucose, simulation_duration=SIMULATION_DURATION):
        #! Not being used, may not need
        """Generator that yields the cell's state after each time step."""
        self.simulation_time = 0
        total_atp_produced = 0
        glucose_processed = 0

        while (
            glucose_processed < glucose and self.simulation_time < simulation_duration
        ):
            if self.metabolites["oxygen"].quantity <= 0:
                self.logger.warning("Oxygen depleted. Stopping simulation.")
                break

            # Calculate ATP at the start of the iteration
            atp_start = (
                self.metabolites["atp"].quantity + self.metabolites["atp"].quantity
            )

            # Check ADP availability and transfer if needed
            if self.metabolites["adp"].quantity < 10:
                self.logger.warning(
                    "Low ADP levels in mitochondrion. Transferring ADP from cytoplasm."
                )
                adp_transfer = min(
                    50, max(0, self.cytoplasm.metabolites["adp"].quantity)
                )  # Transfer up to 50 ADP, but not less than 0
                self.metabolites["adp"].quantity += adp_transfer
                self.metabolites["adp"].quantity -= adp_transfer

            # Implement feedback activation
            adp_activation_factor = 1 + self.cytoplasm.metabolites["adp"].quantity / 500
            self.glycolysis_rate = self.base_glycolysis_rate * adp_activation_factor

            # Glycolysis with updated rate
            glucose_consumed = min(
                1 * self.glycolysis_rate, self.metabolites["glucose"].quantity
            )
            pyruvate = self.glycolysis.perform_glycolysis(glucose_consumed)
            glucose_processed += glucose_consumed

            cytoplasmic_nadh = self.metabolites["nadh"].quantity

            # NADH shuttle
            mitochondrial_nadh = self.mitochondrion.transfer_cytoplasmic_nadh(
                cytoplasmic_nadh
            )

            # Cellular respiration in mitochondrion
            #! Pausing for now
            # mitochondrial_atp = self.mitochondrion.cellular_respiration(pyruvate)

            # Transfer excess ATP from mitochondrion to cytoplasm
            atp_transfer = max(0, self.metabolites["atp"].quantity - 100)
            self.metabolites["atp"].quantity += atp_transfer
            self.metabolites["atp"].quantity -= atp_transfer

            # Calculate ATP at the end of the iteration
            atp_end = (
                self.metabolites["atp"].quantity + self.metabolites["atp"].quantity
            )

            # Calculate ATP produced in this iteration
            delta_atp = atp_end - atp_start
            total_atp_produced += delta_atp

            # Yield the current state of the cell
            yield self.get_cell_state(glucose_processed, total_atp_produced)

            self.simulation_time += self.time_step

        # Final yield after the simulation completes
        yield self.get_cell_state(glucose_processed, total_atp_produced)

    def get_cell_state(self, glucose_processed, total_atp_produced):
        """Helper method to get the current state of the cell."""
        return {
            "simulation_time": self.simulation_time,
            "glucose_processed": glucose_processed,
            "total_atp_produced": total_atp_produced,
            "cytoplasm_atp": self.metabolites["atp"].quantity,
            "mitochondrion_atp": self.metabolites["atp"].quantity,
            "cytoplasm_nadh": self.metabolites["nadh"].quantity,
            "mitochondrion_nadh": self.metabolites["nadh"].quantity,
            "mitochondrion_fadh2": self.metabolites["fadh2"].quantity,
            "mitochondrial_calcium": self.metabolites["calcium"].quantity,
            "proton_gradient": self.mitochondrion.proton_gradient,
            "oxygen_remaining": self.metabolites["oxygen"].quantity,
        }

    def reset(self):
        """Reset the entire cell state."""
        self.atp = 300  # Reset ATP to initial value
        self.metabolites.reset()
        self.simulation_time = 0
        self.logger.info("Cell state reset")

    def check_adenine_nucleotide_balance(self):
        """Check if the total adenine nucleotides have changed."""
        current_total = (
            self.metabolites["ATP"].quantity
            + self.metabolites["ADP"].quantity
            + self.metabolites["AMP"].quantity
        )
        if abs(current_total - self.initial_adenine_nucleotides) > 1e-6:
            self.logger.warning(
                f"Adenine nucleotide imbalance detected. "
                f"Initial: {self.initial_adenine_nucleotides}, "
                f"Current: {current_total}"
            )
        return current_total

    def initialize_simulation(self):
        """Initialize the simulation and record initial adenine nucleotide levels."""
        self.initial_adenine_nucleotides = (
            self.metabolites["ATP"].quantity
            + self.metabolites["ADP"].quantity
            + self.metabolites["AMP"].quantity
        )
        self.logger.info(
            f"Initial total adenine nucleotides: {self.initial_adenine_nucleotides}"
        )
