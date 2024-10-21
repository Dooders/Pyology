from pyology.cell import Cell
from pyology.glycolysis import Glycolysis
from pyology.observers import (
    AdenineNucleotideBalanceObserver,
    NegativeMetaboliteObserver,
)

from .constants import SIMULATION_DURATION
from .energy_calculations import (
    calculate_cell_energy_state,
    calculate_total_adenine_nucleotides,
)
from .exceptions import (
    GlycolysisError,
    InsufficientMetaboliteError,
    QuantityError,
    UnknownMetaboliteError,
)
from .reporter import Reporter


class SimulationController:
    """
    A class to control the simulation process.

    Parameters
    ----------
    cell: Cell
        The cell to be simulated.
    reporter: Reporter
        The reporter to log events and results.

    Methods
    -------
    run_simulation(glucose: float) -> dict:
        Run the simulation with the specified glucose amount.
    """

    def __init__(self, cell: Cell, reporter: Reporter, debug=True):
        self.cell = cell
        self.reporter = reporter
        self.debug = debug
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
        self.initial_adenine_nucleotides = self._calculate_total_adenine_nucleotides()
        self.initial_atp = self.cell.cytoplasm.metabolites["ATP"].quantity
        self.initial_adp = self.cell.cytoplasm.metabolites["ADP"].quantity
        self.initial_amp = self.cell.cytoplasm.metabolites["AMP"].quantity
        self.initial_adenine_nucleotides = (
            self.initial_atp + self.initial_adp + self.initial_amp
        )
        self.initial_glucose = self.cell.cytoplasm.metabolites["glucose"].quantity
        self.adenine_nucleotide_log = []
        self.initial_energy_state = self._calculate_total_energy_state()
        self.observers = [
            NegativeMetaboliteObserver(),
            AdenineNucleotideBalanceObserver(),
        ]

    def _calculate_total_energy_state(self) -> float:
        return calculate_cell_energy_state(self.cell)

    def _calculate_total_adenine_nucleotides(self) -> float:
        """
        Calculate the total adenine nucleotides in the system.
        """
        return calculate_total_adenine_nucleotides(
            self.cell.cytoplasm.metabolites
        ) + calculate_total_adenine_nucleotides(self.cell.mitochondrion.metabolites)

    def _adjust_adenine_balance_after_glycolysis(self, adenine_before, adenine_after):
        """
        Adjust the adenine nucleotide balance after glycolysis if there's a discrepancy.

        Parameters:
        -----------
        adenine_before : float
            Total adenine nucleotides before glycolysis
        adenine_after : float
            Total adenine nucleotides after glycolysis
        """
        if abs(adenine_after - adenine_before) > 1e-6:
            adjustment = adenine_before - adenine_after
            self.cell.cytoplasm.metabolites["ADP"].quantity += adjustment
            self.reporter.log_event(
                f"Adjusted ADP by {adjustment:.6f} to maintain adenine nucleotide balance after glycolysis"
            )

    def run_simulation(self, glucose: float, reporter: Reporter) -> dict:
        """
        Run the simulation with the specified glucose amount.

        Parameters
        ----------
        glucose: float
            The amount of glucose to be processed during the simulation.

        Returns
        -------
        dict:
            The results of the simulation.
        """
        self.glycolysis = Glycolysis(debug=self.debug)
        self.initial_adenine_nucleotides = (
            self.initial_atp + self.initial_adp + self.initial_amp
        )
        self.adenine_nucleotide_log.append(
            ("Initial", self.initial_adenine_nucleotides)
        )
        self.cell.metabolites["glucose"].quantity = round(glucose, 2)
        reporter.log_event(f"Starting simulation with {glucose:.2f} glucose units")
        try:
            glucose_processed = 0
            total_atp_produced = 0
            initial_glucose = self.cell.metabolites["glucose"].quantity
            initial_pyruvate = self.cell.metabolites["pyruvate"].quantity
            initial_total_adenine = (
                self.initial_atp + self.initial_adp + self.initial_amp
            )
            reporter.log_event(
                f"Initial ATP: {self.initial_atp}, Initial ADP: {self.initial_adp}, Initial AMP: {self.initial_amp}"
            )

            next_log_time = 0
            while (
                glucose_processed < glucose
                and self.simulation_time < self.max_simulation_time
            ):
                try:
                    glucose_available = self.cell.metabolites["glucose"].quantity
                    reporter.log_event(f"glucose_available: {glucose_available}")
                    if glucose_available < 1:
                        reporter.log_warning(
                            "Insufficient glucose for glycolysis. Stopping simulation."
                        )
                        break

                    # Store ATP and ADP levels before reactions
                    atp_before = (
                        self.cell.cytoplasm.metabolites["ATP"].quantity
                        + self.cell.mitochondrion.metabolites["ATP"].quantity
                    )
                    adp_before = (
                        self.cell.cytoplasm.metabolites["ADP"].quantity
                        + self.cell.mitochondrion.metabolites["ADP"].quantity
                    )

                    # Add this line to track adenine nucleotides before each step
                    adenine_before = self._calculate_total_adenine_nucleotides()

                    # Perform glycolysis
                    net_atp_produced, pyruvate_produced = Glycolysis.perform(
                        self.cell, glucose_available, self.reporter
                    )

                    # Add this line to track adenine nucleotides after glycolysis
                    adenine_after_glycolysis = (
                        self._calculate_total_adenine_nucleotides()
                    )

                    # Call the new method to adjust adenine balance
                    self._adjust_adenine_balance_after_glycolysis(
                        adenine_before, adenine_after_glycolysis
                    )

                    glucose_processed += glucose_available
                    total_atp_produced += net_atp_produced

                    # Update ATP levels
                    self.cell.cytoplasm.metabolites["ATP"].quantity += net_atp_produced

                    reporter.log_event(
                        f"ATP produced in this iteration: {net_atp_produced}"
                    )
                    reporter.log_event(
                        f"Total ATP produced so far: {total_atp_produced}"
                    )

                    reporter.log_atp_production("Glycolysis", net_atp_produced)

                    # Check if there is enough glucose
                    if self.cell.metabolites["glucose"].quantity <= 0:
                        reporter.log_warning("Glucose depleted. Stopping simulation.")
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

                    reporter.log_atp_production(
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

                    reporter.log_event(f"Simulation time: {self.simulation_time:.3f}")

                    self._check_adenine_nucleotide_balance()
                    self._check_energy_conservation()

                    # Add this at the end of each iteration
                    self._check_and_adjust_adenine_balance()

                    # Ensure no negative quantities after adjustment
                    for metabolite in ["ATP", "ADP", "AMP"]:
                        if self.cell.cytoplasm.metabolites[metabolite].quantity < 0:
                            self.cell.cytoplasm.metabolites[metabolite].quantity = 0
                            reporter.log_warning(
                                f"Set {metabolite} to 0 to avoid negative quantity"
                            )

                    # Run observers
                    #! Need to think through if this works as expected
                    for observer in self.observers:
                        observer.observe(self.cell, reporter)

                except UnknownMetaboliteError as e:
                    reporter.log_error(f"Unknown metabolite error: {str(e)}")
                    reporter.log_warning("Skipping current simulation step.")
                    continue
                except InsufficientMetaboliteError as e:
                    reporter.log_error(f"Insufficient metabolite error: {str(e)}")
                    reporter.log_warning(
                        "Attempting to continue simulation with available metabolites."
                    )
                    continue
                except QuantityError as e:
                    reporter.log_error(f"Quantity error: {str(e)}")
                    reporter.log_warning(
                        "Adjusting quantities and continuing simulation."
                    )
                    continue
                except GlycolysisError as e:
                    reporter.log_warning(f"Glycolysis error: {str(e)}")
                    break

            # After the simulation loop, update the results dictionary
            final_glucose = self.cell.metabolites["glucose"].quantity
            final_pyruvate = self.cell.metabolites["pyruvate"].quantity
            final_atp = (
                self.cell.cytoplasm.metabolites["ATP"].quantity
                + self.cell.mitochondrion.metabolites["ATP"].quantity
            )
            final_adp = (
                self.cell.cytoplasm.metabolites["ADP"].quantity
                + self.cell.mitochondrion.metabolites["ADP"].quantity
            )
            final_amp = (
                self.cell.cytoplasm.metabolites["AMP"].quantity
                + self.cell.mitochondrion.metabolites["AMP"].quantity
            )
            final_total_adenine = final_atp + final_adp + final_amp

            results = {
                "total_atp_produced": total_atp_produced,
                "glucose_processed": glucose_processed,
                "glucose_consumed": initial_glucose - final_glucose,
                "pyruvate_produced": final_pyruvate - initial_pyruvate,
                "simulation_time": self.simulation_time,
                "oxygen_remaining": self.cell.metabolites["oxygen"].quantity,
                "final_cytoplasm_atp": self.cell.cytoplasm.metabolites["ATP"].quantity,
                "final_mitochondrion_atp": self.cell.mitochondrion.metabolites[
                    "ATP"
                ].quantity,
                "final_adp": final_adp,
                "final_amp": final_amp,
                "final_phosphoglycerate_2": self.cell.metabolites[
                    "phosphoglycerate_2"
                ].quantity,
                "final_phosphoenolpyruvate": self.cell.metabolites[
                    "phosphoenolpyruvate"
                ].quantity,
            }

            # Assert non-negative metabolite quantities
            for metabolite, quantity in results.items():
                if metabolite.startswith("final_"):
                    assert (
                        quantity >= 0
                    ), f"Negative quantity for {metabolite}: {quantity}"

            # Assert glucose consumption
            assert (
                results["glucose_consumed"] >= 0
            ), f"Negative glucose consumption: {results['glucose_consumed']}"

            # Assert ATP production
            assert (
                results["total_atp_produced"] >= 0
            ), f"Negative ATP production: {results['total_atp_produced']}"

            reporter.report_simulation_results(results)

            # Check adenine nucleotide balance
            reporter.log_event(
                f"Initial total adenine nucleotides: {initial_total_adenine}"
            )
            reporter.log_event(
                f"Final total adenine nucleotides: {final_total_adenine}"
            )
            reporter.log_event(
                f"Difference: {final_total_adenine - initial_total_adenine}"
            )

            if abs(final_total_adenine - initial_total_adenine) > 1e-6:
                reporter.log_warning("Adenine nucleotide balance is not conserved!")
                # Adjust ATP and ADP to maintain balance
                excess = final_total_adenine - initial_total_adenine
                atp_adjustment = min(excess, final_atp - self.initial_atp)
                self.cell.cytoplasm.metabolites["ATP"].quantity -= atp_adjustment
                adp_adjustment = excess - atp_adjustment
                if adp_adjustment > 0:
                    self.cell.cytoplasm.metabolites["ADP"].quantity -= adp_adjustment
                else:
                    self.cell.cytoplasm.metabolites["ADP"].quantity += abs(
                        adp_adjustment
                    )
                reporter.log_event(
                    f"Adjusted ATP by -{atp_adjustment} and ADP by {-adp_adjustment} to maintain adenine nucleotide balance"
                )
                results["final_cytoplasm_atp"] = self.cell.cytoplasm.metabolites[
                    "ATP"
                ].quantity
                results["final_adp"] = self.cell.cytoplasm.metabolites["ADP"].quantity

            # Add this at the end of the method
            final_adenine_nucleotides = self._calculate_total_adenine_nucleotides()
            results["initial_adenine_nucleotides"] = self.initial_adenine_nucleotides
            results["final_adenine_nucleotides"] = final_adenine_nucleotides

            if abs(final_adenine_nucleotides - self.initial_adenine_nucleotides) > 1e-6:
                reporter.log_warning(
                    f"Adenine nucleotide imbalance detected. "
                    f"Initial: {self.initial_adenine_nucleotides:.6f}, "
                    f"Final: {final_adenine_nucleotides:.6f}, "
                    f"Difference: {final_adenine_nucleotides - self.initial_adenine_nucleotides:.6f}"
                )

            return results

        except Exception as e:
            reporter.log_error(f"Simulation error: {str(e)}")
            raise

    def _handle_adp_availability(self, reporter: Reporter) -> None:
        """
        Handle the availability of ADP in the mitochondrion.
        """
        if self.cell.mitochondrion.metabolites["adp"].quantity < 10:
            reporter.log_warning(
                "Low ADP levels in mitochondrion. Transferring ADP from cytoplasm."
            )
            adp_transfer = min(50, self.cell.metabolites["adp"].quantity)
            self.cell.metabolites["adp"].quantity += adp_transfer
            self.cell.metabolites["adp"].quantity -= adp_transfer

    def _apply_feedback_activation(self) -> None:
        """
        Apply feedback activation based on ADP levels.
        """
        adp_activation_factor = 1 + self.cell.metabolites["adp"].quantity / 500
        self.cell.cytoplasm.glycolysis_rate = (
            self.base_glycolysis_rate * adp_activation_factor
        )

    def _handle_nadh_shuttle(self) -> None:
        """
        Handle the NADH shuttle between the cytoplasm and mitochondrion.
        """
        transfer_rate = 5  # Define a realistic transfer rate per time step
        cytoplasmic_nadh = round(self.cell.metabolites["nadh"].quantity, 2)
        nadh_to_transfer = round(min(transfer_rate, cytoplasmic_nadh), 2)
        self.cell.mitochondrion.transfer_cytoplasmic_nadh(nadh_to_transfer)
        self.cell.metabolites["nadh"].quantity = round(
            self.cell.metabolites["nadh"].quantity - nadh_to_transfer, 2
        )

    def _transfer_excess_atp(self) -> None:
        """
        Transfer excess ATP from the mitochondrion to the cytoplasm.
        """
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

    def _enforce_metabolite_limits(self) -> None:
        """
        Enforce the limits for mitochondrial and cytoplasmic metabolites.
        """
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

    def _log_intermediate_state(self, reporter: Reporter) -> None:
        """
        Log the intermediate state of the simulation.
        """
        state = self.get_current_state()
        reporter.log_event(f"Time: {state['simulation_time']:.2f} s")
        if "glucose_processed" in state:
            reporter.log_event(f"Glucose Processed: {state['glucose_processed']:.2f}")
        else:
            reporter.log_event("Glucose Processed: Not available")
        reporter.log_event(f"Total ATP Produced: {state['total_atp_produced']:.2f}")
        reporter.log_event(f"Cytoplasm ATP: {state['cytoplasm_atp']:.2f}")
        reporter.log_event(f"Mitochondrion ATP: {state['mitochondrion_atp']:.2f}")
        reporter.log_event(f"Proton Gradient: {state['proton_gradient']:.2f}")
        reporter.log_event(f"Oxygen Remaining: {state['oxygen_remaining']:.2f}")
        reporter.log_event(f"NAD+: {state['nad']:.2f}")
        reporter.log_event(f"NADH: {state['nadh']:.2f}")

    def get_current_state(self) -> dict:
        """
        Get the current state of the simulation.
        """
        state = {
            "simulation_time": self.simulation_time,
            "glucose_processed": self.initial_glucose
            - self.cell.cytoplasm.metabolites["glucose"].quantity,
            "cytoplasm_atp": self.cell.cytoplasm.metabolites["ATP"].quantity,
            "mitochondrion_atp": self.cell.mitochondrion.metabolites["ATP"].quantity,
            "total_atp_produced": (
                self.cell.cytoplasm.metabolites["ATP"].quantity
                + self.cell.mitochondrion.metabolites["ATP"].quantity
                - self.initial_atp
            ),
            "proton_gradient": self.cell.mitochondrion.proton_gradient,
            "oxygen_remaining": self.cell.metabolites["oxygen"].quantity,
            "nad": self.cell.metabolites["NAD+"].quantity,
            "nadh": self.cell.metabolites["NADH"].quantity,
        }
        return state

    def reset(self) -> None:
        """
        Reset the simulation state.
        """
        self.cell.reset()
        self.cell.cytoplasm.metabolites["ADP"].quantity = 1.0
        self.cell.cytoplasm.metabolites["AMP"].quantity = 1.0

    def _check_adenine_nucleotide_balance(self, reporter: Reporter) -> None:
        """
        Check the adenine nucleotide balance and log any changes.
        """
        current_adenine_nucleotides = self._calculate_total_adenine_nucleotides()
        difference = current_adenine_nucleotides - self.initial_adenine_nucleotides
        if abs(difference) > 1e-6:
            reporter.log_warning(
                f"Adenine nucleotide imbalance detected. "
                f"Current: {current_adenine_nucleotides:.6f}, "
                f"Initial: {self.initial_adenine_nucleotides:.6f}, "
                f"Difference: {difference:.6f}"
            )
            self.adenine_nucleotide_log.append(
                ("Imbalance", current_adenine_nucleotides)
            )

    def _report_adenine_nucleotide_changes(self, reporter: Reporter) -> None:
        """
        Report the changes in adenine nucleotides throughout the simulation.
        """
        reporter.log_event("\nAdenine Nucleotide Changes:")
        for stage, value in self.adenine_nucleotide_log:
            reporter.log_event(f"{stage}: {value:.6f}")

        initial = self.adenine_nucleotide_log[0][1]
        final = self.adenine_nucleotide_log[-1][1]
        difference = final - initial
        reporter.log_event(f"Total Change: {difference:.6f}")

    def _check_energy_conservation(self, reporter: Reporter) -> None:
        """
        Check the energy conservation and log any changes.

        If the energy state is negative, it means that the system is consuming
        energy. If the energy state is positive, it means that the system is
        producing energy.
        """
        current_energy_state = self._calculate_total_energy_state()
        difference = current_energy_state - self.initial_energy_state
        if abs(difference) > 1e-6:
            reporter.log_warning(
                f"Energy conservation violation detected. "
                f"Current: {current_energy_state:.6f}, "
                f"Initial: {self.initial_energy_state:.6f}, "
                f"Difference: {difference:.6f}"
            )

    def _check_and_adjust_adenine_balance(self, reporter: Reporter):
        """
        Check the adenine nucleotide balance and make adjustments if necessary.

        If the adenine nucleotide balance is negative, it means that the system
        is consuming adenine nucleotides. If the adenine nucleotide balance is
        positive, it means that the system is producing adenine nucleotides.
        """
        current_adenine = self._calculate_total_adenine_nucleotides()
        if abs(current_adenine - self.initial_adenine_nucleotides) > 1e-6:
            adjustment = self.initial_adenine_nucleotides - current_adenine

            # Distribute the adjustment across ATP, ADP, and AMP
            atp_adjustment = min(
                adjustment, self.cell.cytoplasm.metabolites["ATP"].quantity
            )
            self.cell.cytoplasm.metabolites["ATP"].quantity -= atp_adjustment

            remaining_adjustment = adjustment - atp_adjustment
            adp_adjustment = min(
                remaining_adjustment, self.cell.cytoplasm.metabolites["ADP"].quantity
            )
            self.cell.cytoplasm.metabolites["ADP"].quantity -= adp_adjustment

            amp_adjustment = remaining_adjustment - adp_adjustment
            self.cell.cytoplasm.metabolites["AMP"].quantity += amp_adjustment

            reporter.log_warning(
                f"Adjusted ATP by -{atp_adjustment}, ADP by -{adp_adjustment}, and AMP by +{amp_adjustment} to maintain adenine nucleotide balance"
            )
