from lingpy.tests.util import test_data, WithTempDir


class TestWorkflow(WithTempDir):
    def test_Workflow(self):
        from lingpy.basic.workflow import Workflow

        outfile = self.tmp_path('test')
        wf = Workflow(test_data('KSL.qlc'))
        wf.cognate_detection(export='tsv,html', outfile=str(outfile))
        self.assertTrue(self.tmp_path('test.tsv').exists())
        self.assertTrue(self.tmp_path('test.html').exists())
        wf.cognate_detection(cognate_method='lexstat')
