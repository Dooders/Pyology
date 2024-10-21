from typing import Any, Callable, Dict, List

from pydantic import BaseModel, Field


class CommandData(BaseModel):
    """
    Data model that encapsulates the necessary information for executing a command
    on an object, including the object, command, tracked attributes, arguments, and validations.
    """

    obj: Any = Field(
        ..., description="The object on which the command will be executed."
    )
    command: Callable = Field(
        ..., description="The command to be executed on the object."
    )
    tracked_attributes: List[str] = Field(
        ...,
        description="List of attribute names to track before and after the command execution.",
    )
    args: List[Any] = Field(
        default_factory=list,
        description="List of positional arguments to pass to the command.",
    )
    kwargs: Dict[str, Any] = Field(
        default_factory=dict,
        description="Dictionary of keyword arguments to pass to the command.",
    )
    validations: List[Callable[[Any], bool]] = Field(
        default_factory=list,
        description="List of validation functions to run after execution.",
    )

    class Config:
        arbitrary_types_allowed = True
        schema_extra = {
            "example": {
                "obj": "SomeObject()",
                "command": "execute",
                "tracked_attributes": ["status", "result"],
                "args": [1, "test"],
                "kwargs": {"option": True},
                "validations": ["is_valid_result()"],
            }
        }
