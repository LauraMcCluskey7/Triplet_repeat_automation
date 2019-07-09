import unittest
from get_triplet_repeats_3_alleles import *


class test_triplet_repeats(unittest.TestCase):


    def test_get_triplets_table (self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        self.assertEqual(len(triplets_table),7)
        self.assertEqual(len(triplets),12)
        self.assertEqual(triplets_table.iloc[0,3],406.0)        
        self.assertEqual(triplets_table.iloc[4,4],288.9)


    def test_match_control_samples_with_references(self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        self.assertEqual(continue_program, "yes")
        self.assertEqual(len(controls),4)
        self.assertEqual(controls.iloc[0,4],'51')
        self.assertEqual(controls.iloc[3,4],'29')


    def test_find_closest_control_peak_to_sample_peaks(self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        triplets_table=find_closest_control_peak_to_sample_peaks(triplets_table,controls)
        self.assertEqual(triplets_table.iloc[0,7],276)
        self.assertEqual(triplets_table.iloc[0,10],'25')
        self.assertEqual(triplets_table.iloc[3,5],376)
        self.assertEqual(triplets_table.iloc[3,8],'51')
        self.assertEqual(triplets_table.iloc[6,5],323)
        self.assertEqual(triplets_table.iloc[6,8],'45')


    def test_get_number_of_triplet_repeats(self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        triplets_table=find_closest_control_peak_to_sample_peaks(triplets_table,controls)
        triplets_table=get_number_of_triplet_repeats(triplets_table)
        self.assertEqual(triplets_table.iloc[0,12],61)
        self.assertEqual(triplets_table.iloc[0,13],21)
        self.assertEqual(triplets_table.iloc[1,12],54)
        self.assertEqual(triplets_table.iloc[4,13],30)
        self.assertEqual(triplets_table.iloc[6,11],42)


    def test_format_columns(self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        triplets_table=find_closest_control_peak_to_sample_peaks(triplets_table,controls)
        triplets_table=get_number_of_triplet_repeats(triplets_table)
        triplets_table=format_columns(triplets_table, controls, "tester1", "gene1")
        self.assertEqual(triplets_table.iloc[0,2],406)
        self.assertEqual(triplets_table.iloc[0,5],61)
        self.assertEqual(triplets_table.iloc[4,4],30)
        self.assertEqual(triplets_table.iloc[6,1],315)
        self.assertEqual(triplets_table.iloc[6,5],55)


if __name__=='__main__':
    unittest.main()
