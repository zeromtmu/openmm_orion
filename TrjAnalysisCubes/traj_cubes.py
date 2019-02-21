import traceback

from cuberecord import OERecordComputeCube

from Standards import Fields

from floereport import FloeReport, LocalFloeReport

from floe.api import parameter

from orionclient.session import in_orion, OrionSession

from orionclient.types import File

from os import environ


class MDFloeReportCube(OERecordComputeCube):
    version = "0.0.0"
    title = "MDFloeReportCube"
    description = """
    This cube is used to generate an Orion floe report
    """
    classification = [["Floe Reports"]]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    upload = parameter.BooleanParameter(
        'upload',
        default=False,
        help_text="Upload floe report to Amazon S3")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.floe_report_dic = dict()

        if in_orion():
            job_id = environ.get('ORION_JOB_ID')
            self.floe_report = FloeReport.start_report("floe_report", job_id=job_id)
        else:
            self.floe_report = LocalFloeReport.start_report("floe_report")

    def process(self, record, port):

        try:

            if not record.has_value(Fields.title):
                raise ValueError("Missing the title field")

            system_title = record.get_value(Fields.title)

            if not record.has_value(Fields.id):
                raise ValueError("Missing the ID field")

            system_id = record.get_value(Fields.id)

            if not record.has_value(Fields.floe_report):
                raise ValueError("Missing the report field for the system {}".format(system_title + "_" + system_id))

            report_string = record.get_value(Fields.floe_report)

            if not record.has_value(Fields.ligand_name):
                raise ValueError("Missing the ligand name field")

            ligand_name = record.get_value(Fields.ligand_name)

            if not record.has_value(Fields.floe_report_svg_lig_depiction):
                raise ValueError("Missing the ligand  depiction field")

            ligand_svg = record.get_value(Fields.floe_report_svg_lig_depiction)

            if not record.has_value(Fields.floe_report_label):
                floe_report_label = ""
            else:
                floe_report_label = record.get_value(Fields.floe_report_label)

            self.floe_report_dic[system_id] = (report_string, ligand_svg, ligand_name, floe_report_label)

            # Upload Floe Report
            if self.opt['upload']:

                if in_orion():
                    session = OrionSession()

                    file_upload = File.upload(session,
                                              "{}.html".format(system_title),
                                              report_string)

                    session.tag_resource(file_upload, "floe_report")

                    job_id = environ.get('ORION_JOB_ID')

                    if job_id:
                        session.tag_resource(file_upload, "Job {}".format(job_id))

            self.success.emit(record)

        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return

    def end(self):

        try:
            self.opt['Logger'].info("....Generating Floe Report")

            index = self.floe_report.create_page("index", is_index=True)

            index_content = """
            <style>
            .grid { 
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
            grid-gap: 20px;
            align-items: stretch;
            }

            .grid a {
            border: 1px solid #ccc;
            padding: 25px
            }

            .grid svg {
            display: block;  
            max-width: 100%;
            }

            .grid p{
            text-align: center;
            }
            </style>
            <main class="grid">
            """
            # Sort the dictionary keys by using the ligand ID
            for key in sorted(self.floe_report_dic.keys()):

                report_string, ligand_svg, ligand_title, label = self.floe_report_dic[key]

                if len(ligand_title) < 15:
                    page_title = ligand_title
                else:
                    page_title = ligand_title[0:13] + '...'

                page = self.floe_report.create_page(page_title, is_index=False)
                page_link = page.get_link()
                page.set_from_string(report_string)

                index_content += """
                <a href='{}'>
                {}
                <p> {} </p>
                </a>
                """.format(page_link, ligand_svg, label)

            index_content += """
            </main>
            """

            index.set_from_string(index_content)

            self.floe_report.finish_report()

        except Exception as e:
            self.opt['Warning'].warn("It was not possible to generate the floe report: {}".format(str(e)))

        return
