"""Management of locations"""

# ----------------------------- License information --------------------------

# This file is part of the helisol python package.
# Copyright (C) 2023 Olivier Vincent

# The helisol package is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# The helisol package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with the helisol python package.
# If not, see <https://www.gnu.org/licenses/>

import helisol
from copy import copy
import json
from json import JSONDecodeError
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

from .config import CONFIG

LOCATIONS_FOLDER = Path(helisol.__file__).parent.parent / "locations"
LOCATIONS_FILES = "perso.json", "global.json"


@dataclass
class Location:
    name: str
    coords: Tuple
    elevation: int = None

    @classmethod
    def parse(cls, user_input):
        """Generate Location object from user input."""
        try:
            user_input.coords  # Location object or equivalent
        except AttributeError:
            pass
        else:
            return user_input

        if user_input is None:
            return Location.load(CONFIG["default location"])

        elif type(user_input) is str:
            return Location.load(user_input)

        try:
            lat, long = [float(x) for x in user_input]
        except (TypeError, ValueError):
            raise ValueError(f"Unknown location: {user_input}")
        else:
            return Location(name="[Unspecified]", coords=(lat, long))

    @classmethod
    def load(cls, name):
        """Create location object from name.

        Parameters
        ----------
        name : str
            name of location as stored in the JSON database

        Examples
        --------
        >>> location = Location.load('Home')
        """
        for filename in LOCATIONS_FILES:
            file = LOCATIONS_FOLDER / filename
            all_data = cls._from_json(file)
            try:
                data = all_data[name]
                return cls(name=name, **data)
            except KeyError:
                pass
        else:
            raise ValueError(f'Location name "{name}" not found in JSON database')

    def save(self, kind="perso"):
        """Save location data into JSON database.

        Parameters
        ---------
        kind : str
            'perso' (default: non-shared locations), or 'global' (shared)

        Examples
        --------
        >>> location.save('global')
        """
        filename = kind + ".json"
        file = LOCATIONS_FOLDER / filename

        data = self._from_json(file)  # load existing data to add to it

        props = copy(vars(self))
        name = props.pop("name")
        data[name] = props

        self._to_json(file, data)

    def remove(cls, name, kind="perso"):
        """Remove location entry from JSON database"""
        filename = kind + ".json"
        file = LOCATIONS_FOLDER / filename
        data = cls._from_json(file)
        try:
            data.pop(name)
        except KeyError:
            raise ValueError(f"No entry named {name} in {filename} database")
        else:
            cls._to_json(file, data)

    @staticmethod
    def _from_json(file):
        """Load python data (dict or list) from json file"""
        try:
            with open(file, "r", encoding="utf8") as f:
                data = json.load(f)
        except (FileNotFoundError, JSONDecodeError):  # no existing file or data
            data = {}
        return data

    @staticmethod
    def _to_json(file, data):
        """Load python data (dict or list) from json file"""
        with open(file, "w", encoding="utf8") as f:
            json.dump(data, f, indent=4, ensure_ascii=False)
