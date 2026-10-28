"""Test the external-result identity validator without claiming a benchmark."""
from pathlib import Path
import sys, unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
from proseg_comparison import validate_and_compare, run

class ProsegValidationTest(unittest.TestCase):
    def setUp(self):
        self.s=[{'transcript_id':'123','cell_id':'A','feature_name':'G','x_location':'1','y_location':'2'},
                {'transcript_id':'124','cell_id':'UNASSIGNED','feature_name':'G','x_location':'2','y_location':'2'}]
        self.t=[{'transcript_id':'123','gene':'G','observed_x':'1','observed_y':'2','assignment':'0'},
                {'transcript_id':'124','gene':'G','observed_x':'2','observed_y':'2','assignment':'0'}]
        self.c=[{'cell':'0','original_cell_id':'A'}]
    def test_valid(self):
        d=validate_and_compare(self.s,self.t,self.c)
        self.assertEqual(d['same_owner'],1);self.assertEqual(d['changed_owner'],1)
        self.assertEqual(d['count_L1_change'],1)
        self.assertFalse(d['independent_truth_available'])
    def test_excluded_not_silently_joined(self):
        d=validate_and_compare(self.s,self.t[:1],self.c)
        self.assertEqual(d['excluded_source_ids'],['124'])
    def test_duplicates(self):
        with self.assertRaisesRegex(ValueError,'duplicate output'):validate_and_compare(self.s,self.t+self.t,self.c)
    def test_unknown_id(self):
        self.t[0]['transcript_id']='999'
        with self.assertRaisesRegex(ValueError,'unknown output'):validate_and_compare(self.s,self.t,self.c)
    def test_wrong_gene(self):
        self.t[0]['gene']='H'
        with self.assertRaisesRegex(ValueError,'gene mismatch'):validate_and_compare(self.s,self.t,self.c)
    def test_unknown_cell(self):
        self.t[0]['assignment']='99'
        with self.assertRaisesRegex(ValueError,'unknown assigned'):validate_and_compare(self.s,self.t,self.c)
    def test_wrong_coordinates(self):
        self.t[0]['observed_x']='5'
        with self.assertRaisesRegex(ValueError,'coordinate mismatch'):validate_and_compare(self.s,self.t,self.c)
    def test_nonfinite_coordinates(self):
        self.t[0]['observed_x']='nan'
        with self.assertRaisesRegex(ValueError,'coordinate mismatch'):validate_and_compare(self.s,self.t,self.c)
    def test_missing_executable(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'executable unavailable'):run('/nonexistent/proseg',Path(d),Path(d)/'out')

if __name__=='__main__':unittest.main()
