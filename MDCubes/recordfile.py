"""
Copyright (C) 2018 OpenEye Scientific Software
"""
import tempfile

from floe.api.cubes import SourceCube, SinkCube
from floe.api.parameter import StringParameter, FileInputParameter

from orionclient.session import OrionSession, in_orion
from orionclient.types import File

from cuberecord.ports import RecordOutputPort, RawRecordInputPort
from datarecord import OEReadRecords


class RecordFileToRecordConverter(SourceCube):

    title = "Record File to Record Converter"
    description = "Reads a record file and converts to records"
    file = FileInputParameter("file", required=True, title="File to use as input")
    success = RecordOutputPort("success")

    def begin(self):
        if in_orion():
            session = OrionSession()
            if not isinstance(self.args.file, dict):
                file = self.args.file
            else:
                file = self.args.file["file"]
            resource = session.get_resource(File, file)
            self.temp = tempfile.NamedTemporaryFile(suffix=resource.name)
            resource.download_to_file(self.temp.name)
            self.args.file = self.temp.name
        self.ifs = open(self.args.file, "rb")

    def __iter__(self):
        for record in OEReadRecords(self.ifs):
            yield record


class RecordsToRecordFileConverter(SinkCube):

    title = "Record to Record File"
    description = """
A writer that writes a stream of records to a record file (i.e. oedb).
                  """

    file_name = StringParameter(
        "file_name",
        required=True,
        description="Name of the file to create from records"
    )
    intake = RawRecordInputPort("intake")

    def begin(self):
        self.temp = tempfile.NamedTemporaryFile()
        self.ofs = open(self.temp.name, "wb")

    def write(self, record, port):
        self.ofs.write(record)

    def end(self):
        self.ofs.close()
        if in_orion():
            session = OrionSession()
            File.upload(session, self.args.file_name, self.temp.name)
