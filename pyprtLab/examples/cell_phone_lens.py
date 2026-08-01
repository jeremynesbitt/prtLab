# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Chipman cell-phone lens example 9."""

from prtlab import example_cell_phone_lens_ex9_system

if __name__ == "__main__":
    system = example_cell_phone_lens_ex9_system()
    print(system.description)
    print("Surface records:", len(system.surfaces))
