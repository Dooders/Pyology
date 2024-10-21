import logging

from pyology.cell import Cell
from pyology.simulation import Reporter, SimulationController

logging.basicConfig(
    level=logging.DEBUG,  # Change this to DEBUG to see all log messages
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

reporter = Reporter()
cell = Cell(logger=reporter)
sim_controller = SimulationController(cell, reporter)
glucose_amounts = [4]

for glucose in glucose_amounts:
    cell.cytoplasm.metabolites["glucose"].quantity = glucose
    cell.cytoplasm.metabolites["glucose"].initial_quantity = glucose
    reporter.log_event(f"\nSimulating ATP production with {glucose} glucose units:")
    initial_glucose = cell.metabolites["glucose"].quantity
    initial_atp = cell.cytoplasm.metabolites["ATP"].quantity
    initial_adp = cell.metabolites["ADP"].quantity
    initial_amp = cell.metabolites["AMP"].quantity

    # Assert initial metabolite levels are non-negative
    assert initial_glucose >= 0, f"Initial glucose level is negative: {initial_glucose}"
    assert initial_atp >= 0, f"Initial ATP level is negative: {initial_atp}"
    assert initial_adp >= 0, f"Initial ADP level is negative: {initial_adp}"
    assert initial_amp >= 0, f"Initial AMP level is negative: {initial_amp}"

    # Calculate initial total adenine nucleotides
    initial_total_adenine = initial_atp + initial_adp + initial_amp

    reporter.log_event("Initial Metabolite Levels:")
    reporter.log_event(f"ATP: {initial_atp:.2f}")
    reporter.log_event(f"ADP: {initial_adp:.2f}")
    reporter.log_event(f"AMP: {initial_amp:.2f}")

    # Pass the reporter to the run_simulation method
    results = sim_controller.run_simulation(glucose, reporter)

    reporter.log_event("\nAdenine Nucleotide Balance:")
    reporter.log_event(
        f"Initial: {results['initial_adenine_nucleotides']:.6f}, "
        f"Final: {results['final_adenine_nucleotides']:.6f}, "
        f"Difference: {results['final_adenine_nucleotides'] - results['initial_adenine_nucleotides']:.6f}"
    )

    final_atp = results["final_cytoplasm_atp"] + results["final_mitochondrion_atp"]
    final_adp = cell.metabolites["ADP"].quantity
    final_amp = cell.metabolites["AMP"].quantity

    reporter.log_event(
        f"Final: ATP: {final_atp:.2f}, "
        f"ADP: {final_adp:.2f}, "
        f"AMP: {final_amp:.2f}"
    )

    # Assert final metabolite levels are non-negative
    assert final_atp >= 0, f"Final ATP level is negative: {final_atp}"
    assert final_adp >= 0, f"Final ADP level is negative: {final_adp}"
    assert final_amp >= 0, f"Final AMP level is negative: {final_amp}"

    total_initial = initial_atp + initial_adp + initial_amp
    total_final = final_atp + final_adp + final_amp

    reporter.log_event(f"Total Initial Adenine Nucleotides: {total_initial:.2f}")
    reporter.log_event(f"Total Final Adenine Nucleotides: {total_final:.2f}")
    reporter.log_event(f"Difference: {total_final - total_initial:.2f}")

    # Assert conservation of adenine nucleotides with a smaller tolerance
    tolerance = 1e-6  # Decreased tolerance
    assert (
        abs(
            results["final_adenine_nucleotides"]
            - results["initial_adenine_nucleotides"]
        )
        < tolerance
    ), (
        f"Adenine nucleotide conservation violated. "
        f"Initial: {results['initial_adenine_nucleotides']:.6f}, "
        f"Final: {results['final_adenine_nucleotides']:.6f}, "
        f"Difference: {results['final_adenine_nucleotides'] - results['initial_adenine_nucleotides']:.6f}"
    )

    # Assert glucose consumption
    final_glucose = cell.metabolites["glucose"].quantity
    assert final_glucose <= initial_glucose, (
        f"Glucose increased during simulation. "
        f"Initial: {initial_glucose}, Final: {final_glucose}"
    )

    sim_controller.reset()

reporter.log_event("Simulation complete.")
