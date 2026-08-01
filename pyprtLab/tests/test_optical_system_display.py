# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Readable plain-text and notebook display for optical prescriptions."""

from prtlab import quarter_wave_plate_system


def test_optical_system_plain_text_table(capsys):
    system = quarter_wave_plate_system()
    table = system.to_table()
    assert "Wavelength: 0.633 um" in table
    assert "Radius" in table
    assert "uniaxial" in table
    assert "nO: 1.656" in table
    assert "opticAxis" in table

    system.print_table()
    assert capsys.readouterr().out.rstrip() == table.rstrip()


def test_optical_system_html_representation_is_escaped():
    system = quarter_wave_plate_system()
    system.description = "Quarter <wave> plate"
    html = system._repr_html_()
    assert "<table" in html
    assert "Quarter &lt;wave&gt; plate" in html
    assert "<th>Thickness</th>" in html
    assert "opticAxis" in html
