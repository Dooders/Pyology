import logging

from pyology.glycolysis import Glycolysis

from .constants import SIMULATION_DURATION
from .enzymes import Enzyme
from .exceptions import (
    GlycolysisError,
    InsufficientMetaboliteError,
    QuantityError,
    UnknownMetaboliteError,
)
from .reaction import Reaction


class Reporter:
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.atp_production_log = []

    def log_event(self, message):
        self.logger.info(message)

    def log_warning(self, message):
        self.logger.warning(message)

    def log_error(self, message):
        self.logger.error(message)

    def log_atp_production(self, step, atp_produced):
        self.atp_production_log.append((step, atp_produced))
        self.log_event(f"ATP produced in {step}: {atp_produced}")

    def report_simulation_results(self, results):
        self.log_event(
            f"Simulation completed in {results['simulation_time']:.2f} seconds"
        )
        self.log_event(f"Total ATP produced: {results['total_atp_produced']:.2f}")
        self.log_event(f"Glucose processed: {results['glucose_processed']:.2f}")
        self.log_event(f"Oxygen remaining: {results['oxygen_remaining']:.2f}")
        self.log_event(f"Final cytoplasm ATP: {results['final_cytoplasm_atp']:.2f}")
        self.log_event(
            f"Final mitochondrion ATP: {results['final_mitochondrion_atp']:.2f}"
        )

        self.log_event("\nATP Production Breakdown:")
        for step, atp in self.atp_production_log:
            self.log_event(f"  {step}: {atp:.2f}")

        self.atp_production_log.clear()  # Clear the log for the next simulation


class SimulationController:
    def __init__(self, cell, reporter):
        self.cell = cell
        self.reporter = reporter
        self.simulation_duration = SIMULATION_DURATION
        self.simulation_time = 0
        self.time_step = 0.001  # Decreased time step
        self.base_glycolysis_rate = (
            self.cell.cytoplasm.glycolysis_rate
        )  # Store base rate
        self.max_mitochondrial_atp = 100
        self.max_cytoplasmic_atp = 500
        self.max_mitochondrial_nadh = 50
        self.max_cytoplasmic_nadh = 100
        self.max_simulation_time = 20  # Increased max simulation time

    def run_simulation(self, glucose):
        self.cell.initialize_simulation()
        self.initial_adenine_nucleotides = (
            self.cell.metabolites["ATP"].quantity
            + self.cell.metabolites["ADP"].quantity
            + self.cell.metabolites["AMP"].quantity
        )
        self.cell.metabolites["glucose"].quantity = round(glucose, 2)
        self.reporter.log_event(f"Starting simulation with {glucose:.2f} glucose units")
        try:
            glucose_processed = 0
            total_atp_produced = 0
            next_log_time = 10  # Initialize next_log_time here

            while (
                glucose_processed < glucose
                and self.simulation_time < self.max_simulation_time
            ):
                try:
                    glucose_available = self.cell.metabolites["glucose"].quantity
                    self.reporter.log_event(f"glucose_available: {glucose_available}")
                    if glucose_available < 1:
                        self.reporter.log_warning(
                            "Insufficient glucose for glycolysis. Stopping simulation."
                        )
                        break

                    # Store ATP levels before reactions
                    atp_before = (
                        self.cell.cytoplasm.metabolites["atp"].quantity
                        + self.cell.mitochondrion.metabolites["atp"].quantity
                    )

                    # Perform glycolysis
                    pyruvate_produced = Glycolysis.perform(
                        self.cell.cytoplasm, glucose_available
                    )
                    glucose_processed += glucose_available

                    # Calculate ATP produced in this iteration
                    atp_after = (
                        self.cell.cytoplasm.metabolites["atp"].quantity
                        + self.cell.mitochondrion.metabolites["atp"].quantity
                    )
                    atp_produced = round(atp_after - atp_before, 2)
                    total_atp_produced = round(total_atp_produced + atp_produced, 2)

                    self.reporter.log_atp_production("Glycolysis", atp_produced)

                    # Check if there is enough glucose
                    if self.cell.metabolites["glucose"].quantity <= 0:
                        self.reporter.log_warning(
                            "Glucose depleted. Stopping simulation."
                        )
                        break

                    # Check and handle ADP availability
                    self._handle_adp_availability()

                    # Implement feedback activation
                    self._apply_feedback_activation()

                    # Handle NADH shuttle
                    self._handle_nadh_shuttle()

                    # Perform cellular respiration
                    mitochondrial_atp_before = self.cell.mitochondrion.metabolites[
                        "atp"
                    ].quantity
                    #! Pausing for now
                    # mitochondrial_atp = self.cell.mitochondrion.cellular_respiration(pyruvate_produced)

                    mitochondrial_atp_produced = round(
                        self.cell.mitochondrion.metabolites["atp"].quantity
                        - mitochondrial_atp_before,
                        2,
                    )

                    self.reporter.log_atp_production(
                        "Cellular Respiration", mitochondrial_atp_produced
                    )

                    # Transfer excess ATP from mitochondrion to cytoplasm
                    self._transfer_excess_atp()

                    # Ensure metabolite quantities don't exceed limits
                    self._enforce_metabolite_limits()

                    self.simulation_time = round(
                        self.simulation_time + self.time_step, 2
                    )

                    if self.simulation_time >= next_log_time:
                        self._log_intermediate_state()
                        next_log_time = round(
                            next_log_time + 10, 2
                        )  # Schedule next log time

                    self.reporter.log_event(
                        f"Simulation time: {self.simulation_time:.3f}"
                    )

                    self._check_adenine_nucleotide_balance()

                    # Add this at the end of each simulation step
                    self.cell.check_adenine_nucleotide_balance()

                except UnknownMetaboliteError as e:
                    self.reporter.log_error(f"Unknown metabolite error: {str(e)}")
                    self.reporter.log_warning("Skipping current simulation step.")
                    continue
                except InsufficientMetaboliteError as e:
                    self.reporter.log_error(f"Insufficient metabolite error: {str(e)}")
                    self.reporter.log_warning(
                        "Attempting to continue simulation with available metabolites."
                    )
                    continue
                except QuantityError as e:
                    self.reporter.log_error(f"Quantity error: {str(e)}")
                    self.reporter.log_warning(
                        "Adjusting quantities and continuing simulation."
                    )
                    continue
                except GlycolysisError as e:
                    self.reporter.log_warning(f"Glycolysis error: {str(e)}")
                    break
            # print(f"metabolites: {self.cell.metabolites}")
            # Return simulation results
            results = {
                "total_atp_produced": total_atp_produced,
                "glucose_processed": glucose_processed,
                "simulation_time": self.simulation_time,
                "oxygen_remaining": self.cell.metabolites.get("oxygen", 0).quantity,
                "final_cytoplasm_atp": self.cell.metabolites["atp"].quantity,
                "final_mitochondrion_atp": self.cell.metabolites["atp"].quantity,
            }
            self.reporter.report_simulation_results(results)
            return results

        except Exception as e:
            self.reporter.log_error(f"Unhandled simulation error: {str(e)}")
            raise

    def _handle_adp_availability(self):
        if self.cell.mitochondrion.metabolites["adp"].quantity < 10:
            self.reporter.log_warning(
                "Low ADP levels in mitochondrion. Transferring ADP from cytoplasm."
            )
            adp_transfer = min(50, self.cell.metabolites["adp"].quantity)
            self.cell.metabolites["adp"].quantity += adp_transfer
            self.cell.metabolites["adp"].quantity -= adp_transfer

    def _apply_feedback_activation(self):
        adp_activation_factor = 1 + self.cell.metabolites["adp"].quantity / 500
        self.cell.cytoplasm.glycolysis_rate = (
            self.base_glycolysis_rate * adp_activation_factor
        )

    def _handle_nadh_shuttle(self):
        transfer_rate = 5  # Define a realistic transfer rate per time step
        cytoplasmic_nadh = round(self.cell.metabolites["nadh"].quantity, 2)
        nadh_to_transfer = round(min(transfer_rate, cytoplasmic_nadh), 2)
        self.cell.mitochondrion.transfer_cytoplasmic_nadh(nadh_to_transfer)
        self.cell.metabolites["nadh"].quantity = round(
            self.cell.metabolites["nadh"].quantity - nadh_to_transfer, 2
        )

    def _transfer_excess_atp(self):
        atp_excess = max(
            0,
            self.cell.metabolites["atp"].quantity - self.max_mitochondrial_atp,
        )
        transfer_amount = min(
            atp_excess,
            self.max_cytoplasmic_atp - self.cell.metabolites["atp"].quantity,
        )

        self.cell.metabolites["atp"].quantity += transfer_amount
        self.cell.metabolites["atp"].quantity -= transfer_amount

    def _enforce_metabolite_limits(self):
        # Limit mitochondrial metabolites
        self.cell.metabolites["atp"].quantity = min(
            self.cell.metabolites["atp"].quantity,
            self.max_mitochondrial_atp,
        )
        self.cell.metabolites["nadh"].quantity = min(
            self.cell.metabolites["nadh"].quantity,
            self.max_mitochondrial_nadh,
        )

        # Limit cytoplasmic metabolites
        self.cell.metabolites["atp"].quantity = min(
            self.cell.metabolites["atp"].quantity, self.max_cytoplasmic_atp
        )
        self.cell.metabolites["nadh"].quantity = min(
            self.cell.metabolites["nadh"].quantity, self.max_cytoplasmic_nadh
        )

    def _log_intermediate_state(self):
        state = self.get_current_state()
        self.reporter.log_event(f"Time: {state['simulation_time']:.2f} s")
        self.reporter.log_event(f"Glucose Processed: {state['glucose_processed']:.2f}")
        self.reporter.log_event(
            f"Total ATP Produced: {state['total_atp_produced']:.2f}"
        )
        self.reporter.log_event(f"Cytoplasm ATP: {state['cytoplasm_atp']:.2f}")
        self.reporter.log_event(f"Mitochondrion ATP: {state['mitochondrion_atp']:.2f}")
        self.reporter.log_event(f"Proton Gradient: {state['proton_gradient']:.2f}")
        self.reporter.log_event(f"Oxygen Remaining: {state['oxygen_remaining']:.2f}")

    def get_current_state(self):
        return {
            "simulation_time": self.simulation_time,
            "cytoplasm_atp": self.cell.metabolites["atp"].quantity,
            "mitochondrion_atp": self.cell.metabolites["atp"].quantity,
            "cytoplasm_nadh": self.cell.metabolites["nadh"].quantity,
            "mitochondrion_nadh": self.cell.metabolites["nadh"].quantity,
            "mitochondrion_fadh2": self.cell.metabolites["fadh2"].quantity,
            "mitochondrial_calcium": self.cell.metabolites["calcium"].quantity,
            "cytoplasmic_calcium": self.cell.cytoplasmic_calcium.quantity,
            "proton_gradient": self.cell.mitochondrion.proton_gradient,
            "oxygen_remaining": self.cell.metabolites["oxygen"].quantity,
        }

    def reset(self):
        self.cell.reset()  # Assuming Cell class has a reset method
        # Reset any other state variables in SimulationController
        # For example:
        # self.time = 0
        # self.total_atp_produced = 0
        # etc.

    def _check_adenine_nucleotide_balance(self):
        total_adenine_nucleotides = (
            self.cell.metabolites["ATP"].quantity
            + self.cell.metabolites["ADP"].quantity
            + self.cell.metabolites["AMP"].quantity
        )
        if abs(total_adenine_nucleotides - self.initial_adenine_nucleotides) > 1e-6:
            self.reporter.log_warning(
                f"Adenine nucleotide imbalance detected. "
                f"Expected: {self.initial_adenine_nucleotides}, "
                f"Actual: {total_adenine_nucleotides}"
            )
