from threading import Lock
from typing import Dict, List

import yaml

from .exceptions import (
    InsufficientMetaboliteError,
    QuantityError,
    UnknownMetaboliteError,
)

gibbs_free_energies = {
    "ATP": 50,
    "ADP": 30,
    "AMP": 10,
    "GTP": 50,
    "NADH": 158,
    "FADH2": 105,
    "Acetyl_CoA": 31,
    "proton_gradient": 5,
    "glucose": 686,
    "glucose-6-phosphate": 916,
    "fructose-6-phosphate": 916,
    "fructose-1-6-bisphosphate": 1146,
    "glyceraldehyde-3-phosphate": 573,
    "1-3-bisphosphoglycerate": 803,
    "3-phosphoglycerate": 573,
    "2-phosphoglycerate": 573,
    "phosphoenolpyruvate": 803,
    "pyruvate": 343,
    "pyruvate_dehydrogenase": 100,
    "phosphoenolpyruvate_carboxykinase": 100,
    "phosphoenolpyruvate_mutase": 100,
    "pyruvate_kinase": 100,
    "phosphoglycerate_kinase": 100,
    "phosphoglycerate_mutase": 100,
    "phosphoglycerate_phosphatase": 100,
}


class Metabolite:
    """
    A class representing a metabolite in a cell or organelle.

    Attributes
    ----------
    name : str
        The name of the metabolite.
    type : str
        The type of the metabolite.
    quantity : float
        The current quantity of the metabolite.
    max_quantity : float
        The maximum quantity of the metabolite.
    min_quantity : float
        The minimum quantity of the metabolite.
    unit : str
        The unit of the metabolite.
    metadata : dict
        Additional metadata about the metabolite.
    on_change : callable, optional
        A callback function that is called when the quantity of the metabolite changes.
    lock : Lock
        A lock object to ensure thread-safe access to the metabolite's quantity.

    Methods
    -------
    adjust_quantity(amount: float) -> None:
        Adjusts the quantity of the metabolite by the specified amount.
    reset() -> None:
        Resets the quantity of the metabolite to the minimum quantity.
    percentage_filled() -> float:
        Returns the percentage of the metabolite's quantity filled.
    to_dict() -> dict:
        Returns a dictionary representation of the metabolite.
    from_dict(data: dict) -> Metabolite:
        Creates a Metabolite instance from a dictionary.
    __repr__() -> str:
        Returns a string representation of the metabolite.
    """

    def __init__(
        self,
        name: str,
        quantity: float,
        max_quantity: float,
        min_quantity: float = 0,
        unit: str = "mM",
        metadata: dict = None,
        on_change=None,
        type: str = "default",
    ) -> None:
        self.name = name.lower()
        self.label = name
        self.type = type
        self.quantity = float(quantity)
        self.max_quantity = float(max_quantity)
        self.min_quantity = float(min_quantity)
        self.unit = unit
        self.metadata = metadata or {}
        self.on_change = on_change
        self.lock = Lock()

    @property
    def quantity(self):
        return self._quantity

    @quantity.setter
    def quantity(self, value):
        self._quantity = float(value)

    def adjust_quantity(self, amount: float) -> None:
        with self.lock:
            new_quantity = self.quantity + amount
            if new_quantity < self.min_quantity or new_quantity > self.max_quantity:
                raise QuantityError(
                    f"Invalid quantity for {self.name}. Attempted to set to {new_quantity}."
                )
            self.quantity = new_quantity
            if self.on_change:
                self.on_change(self)

    @property
    def energy(self):
        if self.label not in gibbs_free_energies:
            pass
            # raise KeyError(f"Gibbs free energy not found for metabolite: {self.label}")
        return gibbs_free_energies.get(self.label, 1) * self.quantity

    def reset(self) -> None:
        with self.lock:
            self.quantity = self.min_quantity
            if self.on_change:
                self.on_change(self)

    @property
    def percentage_filled(self) -> float:
        return (self.quantity / self.max_quantity) * 100

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "type": self.type,
            "quantity": self.quantity,
            "max_quantity": self.max_quantity,
            "min_quantity": self.min_quantity,
            "unit": self.unit,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: dict):
        return cls(
            data["name"],
            data["quantity"],
            data["max_quantity"],
            data.get("min_quantity", 0),
            data.get("unit", "mM"),
            data.get("metadata", {}),
        )

    def __repr__(self):
        return (
            f"Metabolite(name='{self.name}', quantity={self.quantity}, "
            f"max_quantity={self.max_quantity}, unit='{self.unit}', type='{self.type}')"
        )

    @property
    def state(self) -> dict:
        """
        Returns the current state of the metabolite.

        Returns
        -------
        dict
            A dictionary containing the current state of the metabolite.
        """
        return {
            "name": self.name,
            "label": self.label,
            "type": self.type,
            "quantity": self.quantity,
            "max_quantity": self.max_quantity,
            "min_quantity": self.min_quantity,
            "unit": self.unit,
            "metadata": self.metadata,
        }


class Metabolites:
    """
    A class that manages and stores all metabolites for an organelle.

    This class acts like a dictionary to store Metabolite instances and provides
    methods to manage them collectively.

    Attributes
    ----------
    data : Dict[str, Metabolite]
        A dictionary storing metabolite names as keys and Metabolite instances as values.

    Methods
    -------
    register(name: str, quantity: int, max_quantity: int) -> None:
        Adds a new metabolite or updates an existing one.
    change_quantity(name: str, amount: float) -> None:
        Changes the quantity of a specific metabolite.
    is_available(name: str, amount: float) -> bool:
        Checks if a specific amount of a metabolite is available.
    consume(**metabolites: float) -> None:
        Consumes specified amounts of metabolites.
    produce(**metabolites: float) -> None:
        Produces specified amounts of metabolites.
    validate_all() -> None:
        Validates all metabolites to ensure quantities are within valid ranges.
    """

    DEFAULT_STATE_ATTRIBUTES = ["quantity", "energy"]

    def __init__(self):
        """
        Initializes the Metabolites manager with an empty dictionary.
        """
        self.data: Dict[str, Metabolite] = {}

    def _register(self, name: str, quantity: int, max_quantity: int, metadata: dict = None) -> None:
        """
        Adds a new metabolite or updates an existing one.

        Parameters
        ----------
        name : str
            The name of the metabolite.
        quantity : int
            The initial quantity of the metabolite.
        max_quantity : int
            The maximum allowable quantity of the metabolite.

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
        if name.lower() not in self.data:
            metabolite = Metabolite(name, quantity, max_quantity, metadata=metadata)
            self.data[name.lower()] = metabolite
        else:
            metabolite = self.data[name.lower()]
            new_quantity = min(metabolite.quantity + quantity, metabolite.max_quantity)
            metabolite.quantity = new_quantity

    def register(
        self,
        name: str = None,
        quantity: int = None,
        max_quantity: int = None,
        metadata: dict = None,
        **metabolites,
    ) -> None:
        """
        Adds a new metabolite or updates an existing one. Can accept individual parameters or a dictionary.

        Parameters
        ----------
        name : str, optional
            The name of the metabolite (when registering individually).
        quantity : int, optional
            The initial quantity of the metabolite (when registering individually).
        max_quantity : int, optional
            The maximum allowable quantity of the metabolite (when registering individually).
        metadata : dict, optional
            Additional metadata about the metabolite (when registering individually).
        **metabolites : dict
            A dictionary of metabolites to register, where keys are metabolite names and values are tuples of (quantity, max_quantity).

        Examples
        --------
        # Register individually
        register("glucose", 100, 1000)

        # Register from a dictionary
        register(glucose=(100, 1000), atp=(50, 500))
        """
        if name and quantity is not None and max_quantity is not None:
            self._register(name, quantity, max_quantity, metadata)
        elif metabolites:
            for metabolite_name, (
                metabolite_quantity,
                metabolite_max_quantity,
            ) in metabolites.items():
                self._register(
                    metabolite_name.lower(),
                    metabolite_quantity,
                    metabolite_max_quantity,
                    metadata,
                )
        else:
            raise ValueError(
                "Invalid input. Provide either individual parameters or a dictionary of metabolites."
            )

    def change_quantity(self, name: str, amount: float) -> None:
        """
        Changes the quantity of a specific metabolite.

        Parameters
        ----------
        name : str
            The name of the metabolite.
        amount : float
            The amount to change the quantity by.

        Raises
        ------
        UnknownMetaboliteError
            If the metabolite does not exist.
        QuantityError
            If the new quantity is out of valid range.
        """
        if not isinstance(name, str):
            raise TypeError("Metabolite name must be a string.")
        if not isinstance(amount, (int, float)):
            raise TypeError("Amount must be a number.")
        if name.lower() not in self.data:
            raise UnknownMetaboliteError(f"Unknown metabolite: {name}")

        metabolite = self.data[name.lower()]
        new_quantity = metabolite.quantity + amount

        if new_quantity < metabolite.min_quantity:
            raise QuantityError(
                f"Cannot reduce {name} below {metabolite.min_quantity}. Attempted to set {name} to {new_quantity}."
            )
        if new_quantity > metabolite.max_quantity:
            raise QuantityError(
                f"Cannot exceed max quantity for {name}. Attempted to set {name} to {new_quantity}, but max is {metabolite.max_quantity}."
            )

        metabolite.quantity = new_quantity

    def is_available(self, name: str, amount: float) -> bool:
        """
        Checks if a specific amount of a metabolite is available.

        Parameters
        ----------
        name : str
            The name of the metabolite.
        amount : float
            The amount to check for availability.

        Returns
        -------
        bool
            True if available, False otherwise.

        Raises
        ------
        UnknownMetaboliteError
            If the metabolite does not exist.
        """
        if name.lower() not in self.data:
            raise UnknownMetaboliteError(f"Unknown metabolite: {name}")
        return self.data[name.lower()].quantity >= amount

    def consume(self, **metabolites: float) -> None:
        """
        Consumes specified amounts of metabolites.

        Parameters
        ----------
        metabolites : dict
            Metabolite names and amounts to consume.

        Raises
        ------
        InsufficientMetaboliteError
            If any metabolite is insufficient for consumption.
        """
        # First, check availability
        for name, amount in metabolites.items():
            if not self.is_metabolite_available(name, amount):
                raise InsufficientMetaboliteError(f"Insufficient {name} for reaction")
        # Then, consume them
        for name, amount in metabolites.items():
            self.change_metabolite_quantity(name, -amount)

    def produce(self, **metabolites: float) -> None:
        """
        Produces specified amounts of metabolites.

        Parameters
        ----------
        metabolites : dict
            Metabolite names and amounts to produce.
        """
        for name, amount in metabolites.items():
            if name.lower() not in self.data:
                raise UnknownMetaboliteError(f"Unknown metabolite: {name}")
            self.change_metabolite_quantity(name, amount)

    def validate_all(self) -> None:
        """
        Validates all metabolites to ensure quantities are within valid ranges.

        Raises
        ------
        ValueError
            If any metabolite has an invalid quantity.
        """
        for metabolite in self.data.values():
            if (
                metabolite.quantity < metabolite.min_quantity
                or metabolite.quantity > metabolite.max_quantity
            ):
                raise ValueError(
                    f"Invalid quantity for {metabolite.name}: {metabolite.quantity}"
                )

    @property
    def quantities(self) -> Dict[str, float]:
        return {
            metabolite.name: metabolite.quantity for metabolite in self.data.values()
        }

    def get(self, key: str, default: Metabolite = None) -> Metabolite:
        return self.data.get(key.lower(), default)

    def reset(self):
        for metabolite in self.data.values():
            metabolite.reset()

    def __getitem__(self, key):
        normalized_key = key.lower()
        if normalized_key not in self.data:
            # If the metabolite doesn't exist, create it with default values
            self._register(
                normalized_key, 0, 100
            )  # You may want to adjust these default values
            print(
                f"Warning: Metabolite '{key}' was not found. Created with default values."
            )
        return self.data[normalized_key]

    def __setitem__(self, key: str, value: Metabolite) -> None:
        self.data[key.lower()] = value

    def __delitem__(self, key: str) -> None:
        del self.data[key.lower()]

    def __iter__(self):
        return iter(self.data)

    def __len__(self) -> int:
        return len(self.data)

    def __contains__(self, key: str) -> bool:
        return key.lower() in self.data

    def __repr__(self):
        return f"Metabolites({self.data})"

    def items(self):
        return self.data.items()

    def keys(self):
        return self.data.keys()

    def values(self):
        return self.data.values()

    def __iter__(self):
        return iter(self.data.values())

    @property
    def names(self):
        return [metabolite.name for metabolite in self.data.values()]

    @property
    def labels(self):
        return [metabolite.label for metabolite in self.data.values()]

    @property
    def types(self):
        return [metabolite.type for metabolite in self.data.values()]

    @property
    def quantities(self):
        return [metabolite.quantity for metabolite in self.data.values()]

    @property
    def total_energy(self):
        return sum(metabolite.energy for metabolite in self.data.values())

    @property
    def energy_dict(self):
        return {metabolite.name: metabolite.energy for metabolite in self.data.values()}

    def exists(self, key: str) -> bool:
        """
        Check if a metabolite exists in the collection.

        Parameters
        ----------
        key : str
            The name of the metabolite to check.

        Returns
        -------
        bool
            True if the metabolite exists, False otherwise.
        """
        normalized_key = key.lower()
        return normalized_key in self.data

    def state(self, attributes: list = None) -> dict:
        """
        Returns the current states of all metabolites.

        Parameters
        ----------
        attributes : list, optional
            A list of attribute names to include in the state.
            If None, uses DEFAULT_STATE_ATTRIBUTES.
        no_zero : bool, optional
            If True, only includes metabolites with a quantity greater than zero.

        Returns
        -------
        dict
            A dictionary containing the current states of all metabolites.
        """
        if attributes is None:
            attributes = self.DEFAULT_STATE_ATTRIBUTES

        return {
            x.name: {attr: getattr(x, attr) for attr in attributes if hasattr(x, attr)}
            for x in self.data.values()
        }

    @classmethod
    def _load_metabolite_info(cls, name: str) -> dict:
        """
        Loads the yml file for a metabolite and returns the info as a dict.
        """
        with open(f"pyology/metabolites/{name}.yml", "r", encoding="utf-8") as f:
            return yaml.safe_load(f)

    @classmethod
    def from_list(cls, metabolite_names: List[str]) -> "Metabolites":
        """
        Builds a Metabolites instance from a list of metabolite names.

        Loads the yml file for each metabolite and creates a Metabolite instance for each.
        """
        metabolites = cls()
        if metabolite_names is None:
            metabolite_names = []
        for name in metabolite_names:
            try:
                metabolite_info = cls._load_metabolite_info(name)
                # Extract just the quantity and max_quantity from the YAML data
                quantity = metabolite_info.get("quantity", 0)
                max_quantity = (
                    metabolite_info.get("meta", {})
                    .get("concentration", {})
                    .get("range", {})
                    .get("max", quantity)
                )
                metabolites.register(
                    name,
                    quantity=quantity,
                    max_quantity=max_quantity,
                    metadata=metabolite_info.get("meta", {}),
                )
            except ValueError as e:
                raise ValueError(
                    f"Error processing metabolite '{name}': {str(e)}"
                ) from e
        return metabolites
