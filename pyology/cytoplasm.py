import logging

from utils.command_data import CommandData
from utils.tracking import execute_command

from .glycolysis import Glycolysis
from .organelle import Organelle

logger = logging.getLogger(__name__)


class Cytoplasm(Organelle):
    """
    The cytoplasm is the fluid-filled space inside the cell that contains the
    cell's organelles and performs many of the cell's metabolic functions.

    Attributes
    ----------
    glycolysis_rate (float): The rate at which glucose is processed during glycolysis.

    Methods
    -------
    glycolysis(self, glucose_consumed: float) -> float:
        Perform glycolysis on the specified amount of glucose.
    reset(self) -> None:
        Reset the cytoplasm to its initial state.
    """

    name = "Cytoplasm"

    def __init__(self, logger=None, debug=False):
        super().__init__()
        self.debug = debug
        self.logger = logger or logging.getLogger(__name__)

    def glycolysis(self, glucose_amount: float) -> float:
        #! Not being used, may not need. SHould this me here or somewhere else
        """
        Perform glycolysis on the given amount of glucose.

        Args:
            glucose_amount: The amount of glucose to process.

        Returns:
            The amount of pyruvate produced.
        """
        tracked_attributes = ["glucose", "atp", "adp", "nad", "nadh", "pyruvate"]

        def validate_conservation(obj, initial, final, changes):
            initial_adenine = initial["atp"] + initial["adp"]
            final_adenine = final["atp"] + final["adp"]
            return abs(final_adenine - initial_adenine) < 1e-6

        command_data = CommandData(
            obj=self,
            command=Glycolysis.perform,
            tracked_attributes=tracked_attributes,
            args=[glucose_amount],
            validations=[validate_conservation],
        )

        result = execute_command(command_data, self.logger, self.debug)

        return result["result"]  # This should be the amount of pyruvate produced

    def reset(self) -> None:
        self.__init__()
