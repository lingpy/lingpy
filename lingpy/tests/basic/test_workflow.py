from unittest import TestCase

from lingpy.tests.util import test_data


class TestWorkflow(TestCase):
    def test_Workflow(self):
        from lingpy.basic.workflow import Workflow

        wf = Workflow(test_data('KSL.qlc'))
        wf.cognate_detection()
