import json
import os
from dataclasses import dataclass, field
from typing import Dict, Union, List

from .exceptions import (
    InsufficientMetaboliteError,
    QuantityError,
    UnknownMetaboliteError,
)
from .metabolite import Metabolite, Metabolites


@dataclass
class CellMetabolites:
    """
    A data class that initializes a Metabolites instance with basic cell metabolites.

    Attributes
    ----------
    metabolites : Metabolites
        An instance of the Metabolites class containing all basic metabolites.
    """

    metabolites: Metabolites = field(default_factory=Metabolites)

    def __post_init__(self):
        """
        Post-initialization to add basic metabolites to the Metabolites instance.
        """
        # Load metabolites from JSON file
        json_path = os.path.join(os.path.dirname(__file__), "basic_metabolites.json")
        with open(json_path, "r", encoding="utf-8") as f:
            basic_metabolites = json.load(f)

        # Add each metabolite to the Metabolites instance
        for metabolite_info in basic_metabolites:
            self.metabolites.register(
                name=metabolite_info["name"],
                quantity=metabolite_info["quantity"],
                max_quantity=metabolite_info["max_quantity"],
                type=metabolite_info["type"],
            )


class OrganelleMeta(type):
    """
    Metaclass that enforces required attributes and methods,
    and automatically registers organelle classes.

    Intended to be used as a metaclass for the Organelle class.

    Parameters
    ----------
    name : str
        The name of the organelle class.
    bases : tuple
        The base classes of the organelle class.
    namespace : dict
        The namespace of the organelle class.

    Methods
    -------
    get_registry : classmethod
        Returns the registry of organelle classes.
    """

    _registry = {}  # Class variable to hold registered organelles

    def __new__(cls, name: str, bases: tuple, namespace: dict) -> type:
        """
        Creates a new instance of the class.

        Parameters
        ----------
        name : str
            The name of the organelle class.
        bases : tuple
            The base classes of the organelle class.
        namespace : dict
            The namespace of the organelle class.
        """
        # Check if 'name' is defined in this class or any of its base classes
        if not any("name" in B.__dict__ for B in bases) and "name" not in namespace:
            raise AttributeError(
                f"Class '{name}' must define a 'name' class attribute."
            )

        cls._registry[namespace["name"]] = namespace

        return super().__new__(cls, name, bases, namespace)

    @classmethod
    def get_registry(cls: type) -> dict:
        """
        Returns the registry of organelle classes.

        Returns
        -------
        dict
            The registry of organelle classes.
        """
        return dict(cls._registry)


class Organelle(metaclass=OrganelleMeta):
    """
    Base class for organelles, using OrganelleMeta as the metaclass.

    An organelle is a cellular compartment that contains metabolites and enzymes.

    Attributes
    ----------
    name : str
        The name of the organelle.
    metabolites : dict
        A dictionary of metabolites in the organelle.

    Methods
    -------
    validate_initial_state(self) -> None:
        Validates the initial state of the organelle, including metabolite quantities and rates.
    set_glycolysis_rate(self, rate: float) -> None:
        Sets the glycolysis rate with validation.
    add_metabolite(self, name: str, type: str, quantity: float, max_quantity: float) -> None:
        Adds a metabolite to the organelle or increases its quantity if it already exists.
    change_metabolite_quantity(self, metabolite_name: str, amount: float) -> None:
        Changes the quantity of a metabolite in the organelle.
    is_metabolite_available(self, metabolite: str, amount: float) -> bool:
        Checks if a metabolite is available in the organelle.
    consume_metabolites(self, **metabolites: float) -> None:
        Consumes metabolites from the organelle.
    produce_metabolites(self, **metabolites: float) -> None:
        Produces metabolites in the organelle.
    """

    name = "Organelle"

    def __init__(self):
        self.metabolites = CellMetabolites().metabolites
        self.glycolysis_rate = 1.0
        self.validate_initial_state()

    def validate_initial_state(self) -> None:
        """
        Validates the initial state of the organelle, including metabolite quantities and rates.
        """
        for name, metabolite in self.metabolites.items():
            if metabolite.quantity < 0 or metabolite.quantity > metabolite.max_quantity:
                raise ValueError(
                    f"Invalid initial quantity for {name}: {metabolite.quantity}"
                )

        if self.glycolysis_rate <= 0:
            raise ValueError(f"Invalid glycolysis rate: {self.glycolysis_rate}")

    def set_glycolysis_rate(self, rate: float) -> None:
        """
        Sets the glycolysis rate with validation.

        Parameters
        ----------
        rate : float
            The new glycolysis rate.

        Raises
        ------
        ValueError
            If the rate is not positive.
        """
        if rate <= 0:
            raise ValueError(f"Glycolysis rate must be positive. Got: {rate}")
        self.glycolysis_rate = rate

    def add_metabolite(
        self, name: str, type: str, quantity: float, max_quantity: float
    ) -> None:
        """
        Add a metabolite to the organelle or increase its quantity if it already exists.

        Parameters
        ----------
        name : str
            The name of the metabolite.
        type : str
            The type of the metabolite.
        quantity : float
            The quantity of the metabolite.
        max_quantity : float
            The maximum quantity of the metabolite.

        Raises
        ------
        ValueError
            If the quantity is negative or exceeds the maximum quantity.
        """
        if quantity < 0:
            raise ValueError(f"Quantity must be non-negative. Got: {quantity}")

        if quantity > max_quantity:
            raise ValueError(
                f"Initial quantity {quantity} exceeds max quantity {max_quantity}."
            )

        if name in self.metabolites:
            new_quantity = min(self.metabolites[name].quantity + quantity, max_quantity)
            self.metabolites[name].quantity = new_quantity
        else:
            metabolite = Metabolite(name, quantity, max_quantity)
            metabolite.type = type  # Set the type after creation
            self.metabolites[name] = metabolite

    def change_metabolite_quantity(self, metabolite_name: str, amount: float) -> None:
        """
        Changes the quantity of a metabolite in the organelle.

        Parameters
        ----------
        metabolite_name : str
            The name of the metabolite.
        amount : float
            The amount to change the quantity by.
        """
        if not isinstance(metabolite_name, str):
            raise TypeError("Metabolite name must be a string.")
        if not isinstance(amount, (int, float)):
            raise TypeError("Amount must be a number.")
        if metabolite_name not in self.metabolites:
            raise UnknownMetaboliteError(f"Unknown metabolite: {metabolite_name}")

        metabolite = self.metabolites[metabolite_name]
        new_quantity = metabolite.quantity + amount

        if new_quantity < 0:
            raise QuantityError(
                f"Cannot reduce {metabolite_name} below zero. Attempted to set {metabolite_name} to {new_quantity}."
            )
        if new_quantity > metabolite.max_quantity:
            raise QuantityError(
                f"Cannot exceed max quantity for {metabolite_name}. Attempted to set {metabolite_name} to {new_quantity}, but max is {metabolite.max_quantity}."
            )

        metabolite.quantity = new_quantity

    def is_metabolite_available(self, metabolite: str, amount: float) -> bool:
        """
        Checks if a metabolite is available in the organelle.

        Parameters
        ----------
        metabolite : str
            The name of the metabolite.
        amount : float
            The amount of the metabolite to check for.

        Returns
        -------
        bool
            True if the metabolite is available, False otherwise.
        """
        if metabolite not in self.metabolites:
            raise UnknownMetaboliteError(f"Unknown metabolite: {metabolite}")
        return self.metabolites[metabolite].quantity >= amount

    def consume_metabolites(self, **metabolites: float) -> None:
        """
        Consumes metabolites from the organelle.

        Parameters
        ----------
        metabolites : dict
            A dictionary of metabolites to consume.
        """
        for metabolite, amount in metabolites.items():
            if not self.is_metabolite_available(metabolite, amount):
                raise InsufficientMetaboliteError(
                    f"Insufficient {metabolite} for reaction"
                )
            self.change_metabolite_quantity(metabolite, -amount)

    def produce_metabolites(self, **metabolites: float) -> None:
        """
        Produces metabolites in the organelle.

        Parameters
        ----------
        metabolites : dict
            A dictionary of metabolites to produce.
        """
        for metabolite, amount in metabolites.items():
            self.change_metabolite_quantity(metabolite, amount)

    def get_metabolite_quantity(self, metabolite: Union[str, List[str]]) -> Union[float, List[float]]:
        """
        Returns the quantity of one or more metabolites in the organelle.

        Parameters
        ----------
        metabolite : str or List[str]
            The name of the metabolite or a list of metabolite names.

        Returns
        -------
        float or List[float]
            The quantity of the metabolite(s).

        Raises
        ------
        UnknownMetaboliteError
            If any of the requested metabolites are not found in the organelle.
        """
        if isinstance(metabolite, str):
            if metabolite not in self.metabolites:
                raise UnknownMetaboliteError(f"Unknown metabolite: {metabolite}")
            return self.metabolites[metabolite].quantity
        elif isinstance(metabolite, list):
            quantities = []
            for met in metabolite:
                if met not in self.metabolites:
                    raise UnknownMetaboliteError(f"Unknown metabolite: {met}")
                quantities.append(self.metabolites[met].quantity)
            return quantities
        else:
            raise TypeError("Metabolite must be a string or a list of strings")

    def set_metabolite_quantity(self, metabolite: str, quantity: float) -> None:
        """
        Sets the quantity of a metabolite in the organelle.

        Parameters
        ----------
        metabolite : str
            The name of the metabolite.
        quantity : float
            The quantity of the metabolite.
        """
        if metabolite not in self.metabolites:
            raise UnknownMetaboliteError(f"Unknown metabolite: {metabolite}")
        self.metabolites[metabolite].quantity = quantity

    def get_metabolite(self, metabolite_name):
        """
        Get a metabolite from the cytoplasm.

        Parameters:
        -----------
        metabolite_name : str
            The name of the metabolite to retrieve.

        Returns:
        --------
        Metabolite
            The requested metabolite object.

        Raises:
        -------
        KeyError
            If the metabolite is not found in the cytoplasm.
        """
        return self.metabolites[metabolite_name]
