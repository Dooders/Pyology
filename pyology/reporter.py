import logging
import sys


class Reporter:
    """
    A class to report events and log messages during the simulation.

    Methods
    -------
    log_event(message: str) -> None:
        Log an event message.
    log_warning(message: str) -> None:
        Log a warning message.
    log_error(message: str) -> None:
        Log an error message.
    log_atp_production(step: str, atp_produced: float) -> None:
        Log the ATP production for a specific step.
    report_simulation_results(results: dict) -> None:
        Report the simulation results.
    """

    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        
        # Create console handler and set level to INFO
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        
        # Create formatter
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        # Add formatter to console handler
        console_handler.setFormatter(formatter)
        
        # Add console handler to logger
        self.logger.addHandler(console_handler)
        
        self.atp_production_log = []

    def info(self, message: str) -> None:
        """
        Log an info message.
        """
        self.logger.info(message)

    def log_event(self, message: str) -> None:
        """
        Log an event message.
        """
        self.logger.info(message)

    def log_warning(self, message: str) -> None:
        """
        Log a warning message.
        """
        self.logger.warning(message)

    def log_error(self, message: str) -> None:
        """
        Log an error message.
        """
        self.logger.error(message)

    def log_atp_production(self, step: str, atp_produced: float) -> None:
        """
        Log the ATP production for a specific step.
        """
        self.atp_production_log.append((step, atp_produced))
        self.log_event(f"ATP produced in {step}: {atp_produced}")

    def report_simulation_results(self, results: dict) -> None:
        """
        Report the simulation results.
        """
        self.log_event(
            f"Simulation completed in {results['simulation_time']:.2f} seconds"
        )
        self.log_event(f"Total ATP produced: {results['total_atp_produced']:.2f}")
        self.log_event(f"Glucose processed: {results['glucose_processed']:.2f}")
        self.log_event(f"Glucose consumed: {results['glucose_consumed']:.2f}")
        self.log_event(f"Pyruvate produced: {results['pyruvate_produced']:.2f}")
        self.log_event(f"Oxygen remaining: {results['oxygen_remaining']:.2f}")
        self.log_event(f"Final cytoplasm ATP: {results['final_cytoplasm_atp']:.2f}")
        self.log_event(
            f"Final mitochondrion ATP: {results['final_mitochondrion_atp']:.2f}"
        )
        self.log_event(
            f"2-Phosphoglycerate remaining: {results['final_phosphoglycerate_2']:.2f}"
        )
        self.log_event(
            f"Phosphoenolpyruvate produced: {results['final_phosphoenolpyruvate']:.2f}"
        )

        self.log_event("\nATP Production Breakdown:")
        for step, atp in self.atp_production_log:
            self.log_event(f"  {step}: {atp:.2f}")

        self.atp_production_log.clear()  # Clear the log for the next simulation

    # Add this new method
    def error(self, message: str) -> None:
        """
        Log an error message (alias for log_error).
        """
        self.log_error(message)
