import logging

from pyology.cell import Cell
from pyology.glycolysis import Glycolysis
from pyology.reporter import Reporter


class GlycolysisSimulation:
    def __init__(self, cell: "Cell", debug=True):
        self.cell = cell
        self.debug = debug

    def run(self, glucose_units: float, logger: logging.Logger):
        logger.info("Starting glycolysis simulation")
        # Investment phase
        Glycolysis.investment_phase(self.cell, glucose_units, logger=logger)

        # Yield phase
        Glycolysis.yield_phase(self.cell, glucose_units * 2, logger=logger)


reporter = Reporter()
cell = Cell(logger=reporter)
sim_controller = GlycolysisSimulation(cell, debug=True)
sim_controller.run(4, logger=reporter)
