import logging
from logging import Logger
from typing import Any, Dict

from .command_data import CommandData

default_logger = logging.getLogger(__name__)


def execute_command(
    command_data: CommandData, logger: Logger = default_logger, debug: bool = False
) -> Dict[str, Any]:
    """
    Executes a command on an object based on the provided CommandData object,
    tracks specified attributes, logs the results, and performs optional validations.

    Parameters
    ----------
    command_data : CommandData
        The CommandData object containing all necessary information for command execution.
    logger : Logger, optional
        The logger object to use for logging. Defaults to the default logger.
    debug : bool, optional
        Whether to enable debug logging (default is False).

    Returns
    -------
    Dict[str, Any]
        A dictionary containing initial and final values of tracked attributes,
        their changes, and validation results.
    """
    obj = command_data.obj
    command = command_data.command
    tracked_attributes = command_data.tracked_attributes
    args = command_data.args
    kwargs = command_data.kwargs
    validations = command_data.validations

    # Log initial values of tracked attributes
    initial_values = {}
    for attr in tracked_attributes:
        if hasattr(obj, attr) and not callable(getattr(obj, attr)):
            try:
                value = getattr(obj, attr)
                initial_values[attr] = value
                if debug:
                    logger.debug(f"Initial {attr}: {value}")
            except AttributeError:
                logger.warning(f"Attribute '{attr}' not found on object. Skipping.")

    # Execute the command
    try:
        if callable(command):
            result = command(obj, *args, **kwargs)
        else:
            result = getattr(obj, command)(*args, **kwargs)
        if debug:
            logger.debug(f"Executed command '{command}' with result: {result}")
    except Exception as e:
        logger.error(f"Error executing command '{command}': {str(e)}")
        raise

    # Log and store final values, calculate changes
    final_values = {}
    changes = {}
    for attr in tracked_attributes:
        if attr in initial_values:
            try:
                final_value = getattr(obj, attr)
                final_values[attr] = final_value
                if debug:
                    logger.debug(f"Final {attr}: {final_value}")

                if isinstance(initial_values[attr], (int, float)) and isinstance(
                    final_value, (int, float)
                ):
                    change = final_value - initial_values[attr]
                    changes[attr] = change
                    if debug:
                        logger.debug(f"Change in {attr}: {change}")
                else:
                    if debug:
                        logger.debug(
                            f"{attr} changed from {initial_values[attr]} to {final_value}"
                        )
            except AttributeError:
                logger.warning(
                    f"Attribute '{attr}' not found on object after execution. Skipping."
                )

    # Perform validations
    validation_results = {}
    for idx, validation in enumerate(validations, start=1):
        try:
            validation_result = validation(obj, initial_values, final_values, changes)
            validation_results[f"validation_{idx}"] = validation_result
            if not validation_result:
                logger.warning(f"Validation {idx} failed: {validation.__name__}")
        except Exception as e:
            logger.error(f"Error during validation {validation.__name__}: {str(e)}")
            validation_results[f"validation_{idx}"] = f"error: {str(e)}"

    # Prepare and return results
    return {
        "result": result,
        "initial_values": initial_values,
        "final_values": final_values,
        "changes": changes,
        "validation_results": validation_results,
    }
