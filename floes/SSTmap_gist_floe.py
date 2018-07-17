#!/usr/bin/env python

# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.


from floe.api import WorkFloe
from cuberecord import DataSetWriterCube, DataSetReaderCube
from TrjAnalysisCubes.sstmap import SSTMapGist

job = WorkFloe("SSTMAP GIST")

job.description = """
Testing Floe

Ex. python floes/up.py --ligands ligands.oeb
--ofs-data_out prep.oeb

Parameters:
-----------
ligands (file): OEB file of the prepared ligands

Outputs:
--------
ofs: Output file
"""

job.classification = [
    ["OpenEye", "Example"],
    ["OpenEye", "Custom"]
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DataSetReaderCube("ifs")
ofs = DataSetWriterCube("ofs")
fail = DataSetWriterCube("fail")
ligand_ifs = DataSetReaderCube("ligand_ifs")

# Promotes the parameter
ifs.promote_parameter("data_in", promoted_name="in")
ligand_ifs.promote_parameter("data_in", promoted_name="ligand")

ofs.promote_parameter("data_out", promoted_name="out")
fail.set_parameters(data_out="fail.oedb")

sstmap_gist = SSTMapGist("sstmap_gist", title="SSTMap GIST")

job.add_cubes(ifs, ligand_ifs, ofs, fail, sstmap_gist)

ifs.success.connect(sstmap_gist.intake)
ligand_ifs.success.connect(sstmap_gist.ligand_port)
sstmap_gist.success.connect(ofs.intake)
sstmap_gist.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
