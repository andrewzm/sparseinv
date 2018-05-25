# sparseinv: An R Software package for computing the sparse inverse subset with the Takahashi equations with large datasets.
# Copyright (c) 2017 Andrew Zammit-Mangion
# Author: Andrew Zammit-Mangion, azm (at) uow.edu.au
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

## Unlooad compiled code when unloading FRK
.onUnload <- function (libpath) {
  library.dynam.unload("sparseinv", libpath)
}
