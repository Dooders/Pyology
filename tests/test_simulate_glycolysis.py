# test_simulate_glycolysis.py

from unittest.mock import Mock, patch

import pytest

# Import the GlycolysisSimulation class from your simulation code
from simulate_glycolysis import GlycolysisSimulation


# Create mock versions of the pyology modules
class MockCell:
    def __init__(self, logger):
        self.logger = logger
        self.atp = 10  # Assume starting with 10 ATP molecules


class MockGlycolysis:
    @staticmethod
    def investment_phase(cell, glucose_units, logger=None):
        total_atp_consumed = glucose_units * 2  # Consume 2 ATP per glucose
        cell.atp -= total_atp_consumed
        logger.info(f"ATP consumed in investment phase: {total_atp_consumed}")

    @staticmethod
    def yield_phase(cell, g3p_units, logger=None):
        total_atp_produced = g3p_units * 2  # Produce 2 ATP per G3P
        cell.atp += total_atp_produced
        logger.info(f"ATP produced in yield phase: {total_atp_produced}")


class MockReporter:
    def __init__(self):
        self.logs = []

    def info(self, message):
        self.logs.append(message)


@pytest.fixture
def mock_cell():
    reporter = MockReporter()
    cell = MockCell(logger=reporter)
    return cell


@pytest.fixture
def mock_glycolysis():
    return MockGlycolysis


@pytest.fixture
def mock_reporter():
    return MockReporter()


def test_full_glycolysis_simulation(mock_cell, mock_glycolysis, mock_reporter):
    # Patch the Glycolysis class in simulate_glycolysis module
    with patch("simulate_glycolysis.Glycolysis", mock_glycolysis):
        # Create an instance of GlycolysisSimulation
        sim = GlycolysisSimulation(mock_cell, debug=True)
        # Run the simulation with 4 glucose units
        sim.run(4, logger=mock_reporter)
        # Calculate expected ATP levels
        initial_atp = 10
        atp_consumed = 4 * 2  # Investment phase ATP consumption
        atp_after_investment = initial_atp - atp_consumed
        atp_produced = 8 * 2  # Yield phase ATP production (4 glucose * 2 G3P * 2 ATP)
        final_atp = atp_after_investment + atp_produced
        # Assertions
        assert (
            mock_cell.atp == final_atp
        ), f"Expected final ATP to be {final_atp}, got {mock_cell.atp}"
        assert f"ATP consumed in investment phase: {atp_consumed}" in mock_reporter.logs
        assert f"ATP produced in yield phase: {atp_produced}" in mock_reporter.logs


def test_investment_phase(mock_cell, mock_glycolysis, mock_reporter):
    with patch("simulate_glycolysis.Glycolysis", mock_glycolysis):
        # Run only the investment phase
        mock_glycolysis.investment_phase(mock_cell, 4, logger=mock_reporter)
        # Expected ATP after investment phase
        expected_atp = 10 - (4 * 2)
        assert (
            mock_cell.atp == expected_atp
        ), f"Expected ATP after investment phase to be {expected_atp}, got {mock_cell.atp}"
        assert "ATP consumed in investment phase: 8" in mock_reporter.logs


def test_yield_phase(mock_cell, mock_glycolysis, mock_reporter):
    # Set cell ATP to the value after investment phase
    mock_cell.atp = 2  # Starting ATP after investment phase
    with patch("simulate_glycolysis.Glycolysis", mock_glycolysis):
        # Run the yield phase
        mock_glycolysis.yield_phase(mock_cell, 8, logger=mock_reporter)
        # Expected ATP after yield phase
        expected_atp = 2 + (8 * 2)
        assert (
            mock_cell.atp == expected_atp
        ), f"Expected ATP after yield phase to be {expected_atp}, got {mock_cell.atp}"
        assert "ATP produced in yield phase: 16" in mock_reporter.logs


def test_atp_levels_during_simulation(mock_cell, mock_glycolysis, mock_reporter):
    # Mock methods to track ATP levels after each step
    atp_levels = []

    def investment_phase(cell, glucose_units, logger=None):
        for i in range(int(glucose_units)):
            cell.atp -= 2  # Consume 2 ATP per glucose
            atp_levels.append(cell.atp)
            logger.info(f"ATP after processing glucose unit {i+1}: {cell.atp}")

    def yield_phase(cell, g3p_units, logger=None):
        for i in range(int(g3p_units)):
            cell.atp += 2  # Produce 2 ATP per G3P
            atp_levels.append(cell.atp)
            logger.info(f"ATP after processing G3P unit {i+1}: {cell.atp}")

    # Patch Glycolysis methods
    with patch.object(
        MockGlycolysis, "investment_phase", investment_phase
    ), patch.object(MockGlycolysis, "yield_phase", yield_phase), patch(
        "simulate_glycolysis.Glycolysis", MockGlycolysis
    ):
        sim = GlycolysisSimulation(mock_cell, debug=True)
        sim.run(4, logger=mock_reporter)

    # Expected ATP levels after each glucose unit
    expected_investment_levels = [8, 6, 4, 2]
    # Expected ATP levels after each G3P unit
    expected_yield_levels = [4, 6, 8, 10, 12, 14, 16, 18]
    expected_atp_levels = expected_investment_levels + expected_yield_levels

    assert (
        atp_levels == expected_atp_levels
    ), f"Expected ATP levels to be {expected_atp_levels}, got {atp_levels}"
