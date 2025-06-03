# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

from pathlib import Path


def get_some_list(field):
    """Get a list, whether it's in a file or directly specified."""
    if isinstance(field, str):
        if Path(field).is_file():
            with Path(field).open() as f:
                slist = [line.rstrip() for line in f.readlines()]
        else:
            slist = [field]
    elif isinstance(field, list):
        slist = field

    return slist
