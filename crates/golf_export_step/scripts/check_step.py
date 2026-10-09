"""Checks a STEP file with OCCT, through FreeCAD.

    FreeCADCmd crates/golf_export_step/scripts/check_step.py path/to/file.step

Prints each part's validity, solids, shells, faces and volume, then the result
of OCCT's strict (BOP) check on the whole file. FreeCADCmd passes its own
arguments through, so the file is taken as the last one.
"""

import sys

import FreeCAD
import Import
import Part

path = sys.argv[-1]
doc = FreeCAD.newDocument("check")
Import.insert(path, doc.Name)
for o in doc.Objects:
    shape = getattr(o, "Shape", None)
    if shape is None or shape.isNull() or not shape.Solids:
        continue
    closed = all(shell.isClosed() for shell in shape.Shells)
    print(
        f"{o.Label:<16} valid={shape.isValid()} closed={closed} solids={len(shape.Solids)}"
        f" shells={len(shape.Shells)} faces={len(shape.Faces)} volume={shape.Volume:.3f}"
    )

whole = Part.Shape()
whole.read(path)
try:
    whole.check(True)
    print("strict check: no errors")
except Exception as error:
    print(f"strict check: {error}")
