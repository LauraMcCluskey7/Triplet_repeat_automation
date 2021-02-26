import unittest
from triplet_repeat_automation import *


class test_triplet_repeats(unittest.TestCase):


    def test_get_triplets_table (self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        self.assertEqual(len(triplets_table),7)
        self.assertEqual(len(triplets),12)
        Size1=list(triplets_table["Size 1"])
        Size2=list(triplets_table["Size 2"])
        Size3=list(triplets_table["Size 3"])
        self.assertEqual(Size1, [306, 300, 342, 362, 288, 306, 315]) 
        self.assertEqual(Size2[0], 406.0)
        self.assertEqual(Size2[1], 385.0)
        self.assertEqual(Size2[4], 486.0)
        self.assertEqual(Size2[6], 387.0)
        self.assertEqual(Size3[0], 265.0)
        self.assertEqual(Size3[4], 289.0)



    def test_match_control_samples_with_references(self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        self.assertEqual(continue_program, "yes")
        self.assertEqual(len(controls),4)
        controls_size1=list(controls["Size 1"])
        controls_size2=list(controls["Size 2"])
        controls_triplets_1=list(controls["triplets_1"])
        controls_triplets_2=list(controls["triplets_2"])
        self.assertEqual(controls_size1, [288.7, 297.4, 296.5, 304.5])
        self.assertEqual(controls_size2, [375.8, 323.0, 275.9, 284.1])
        self.assertEqual(controls_triplets_1, ['30', '33', '32', '35'])
        self.assertEqual(controls_triplets_2, ['51', '45', '25', '29'])


        triplets,triplets_table=get_triplets_table("gene1","tester2")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        self.assertEqual(continue_program, "yes")


        triplets,triplets_table=get_triplets_table("gene1","tester3")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        self.assertEqual(continue_program, "no")


    def test_get_closest_value(self):
        #if two values in the array are both the same distance away from x, it will take the first value.
        self.assertEqual(get_closest_value(5, numpy.array([4,9,2,6])),4)
        self.assertEqual(get_closest_value(5, numpy.array([6,9,2,4])),6)
        self.assertEqual(get_closest_value(32,numpy.array([66.5,29.5,23.6,44.4, 34.6])),29.5)


    def test_find_closest_control_peak_to_sample_peaks(self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        triplets_table=find_closest_control_peak_to_sample_peaks(triplets_table,controls)
        size1=list(triplets_table["Size 1"])
        size2=list(triplets_table["Size 2"])
        size3=list(triplets_table["Size 3"])
        closest_1=list(triplets_table["closest_1"])
        closest_2=list(triplets_table["closest_2"])
        closest_3=list(triplets_table["closest_3"])
        repeats_closest_1=list(triplets_table["repeats_closest_1"])
        repeats_closest_2=list(triplets_table["repeats_closest_2"])
        repeats_closest_3=list(triplets_table["repeats_closest_3"])
        self.assertEqual(size1, [306, 300, 342, 362, 288, 306,315])
        self.assertEqual(size2[0], 406.0)
        self.assertEqual(size2[1], 385.0)
        self.assertEqual(size2[4], 486.0)
        self.assertEqual(size2[6], 387.0)
        self.assertEqual(size3[0], 265.0)
        self.assertEqual(size3[4], 289.0)
        self.assertEqual(closest_1, [305, 297, 323, 376, 289, 305, 323])
        self.assertEqual(closest_2[0], 376.0)
        self.assertEqual(closest_2[1], 376.0)
        self.assertEqual(closest_2[4], 376.0)
        self.assertEqual(closest_2[6], 376.0)

        self.assertEqual(closest_3[0], 276.0)
        self.assertEqual(closest_3[4], 289.0)

        self.assertEqual(repeats_closest_1, ['35', '32', '45', '51', '30', '35', '45'])
        self.assertEqual(repeats_closest_2[0], '51')
        self.assertEqual(repeats_closest_2[1], '51')
        self.assertEqual(repeats_closest_2[4], '51')
        self.assertEqual(repeats_closest_2[6], '51')

        self.assertEqual(repeats_closest_3[0], '25')
        self.assertEqual(repeats_closest_3[4], '30')


    def test_get_number_of_triplet_repeats(self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        triplets_table=find_closest_control_peak_to_sample_peaks(triplets_table,controls)
        triplets_table=get_number_of_triplet_repeats(triplets_table)


        Sample=list(triplets_table["Sample"])
        size1=list(triplets_table["Size 1"])
        size2=list(triplets_table["Size 2"])
        size3=list(triplets_table["Size 3"])
        closest_1=list(triplets_table["closest_1"])
        closest_2=list(triplets_table["closest_2"])
        closest_3=list(triplets_table["closest_3"])
        repeats_closest_1=list(triplets_table["repeats_closest_1"])
        repeats_closest_2=list(triplets_table["repeats_closest_2"])
        repeats_closest_3=list(triplets_table["repeats_closest_3"])
        repeats_1=list(triplets_table["Repeats_1"])
        repeats_2=list(triplets_table["Repeats_2"])
        repeats_3=list(triplets_table["Repeats_3"])


        self.assertEqual(Sample, ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6","Sample7"])
        self.assertEqual(size1, [306, 300, 342, 362, 288, 306,315])
        self.assertEqual(size2[0], 406.0)
        self.assertEqual(size2[1], 385.0)
        self.assertEqual(size2[4], 486.0)
        self.assertEqual(size2[6], 387.0)
        self.assertEqual(size3[0], 265.0)
        self.assertEqual(size3[4], 289.0)
        self.assertEqual(closest_1, [305, 297, 323, 376, 289, 305, 323])
        self.assertEqual(closest_2[0], 376.0)
        self.assertEqual(closest_2[1], 376.0)
        self.assertEqual(closest_2[4], 376.0)
        self.assertEqual(closest_2[6], 376.0)

        self.assertEqual(closest_3[0], 276.0)
        self.assertEqual(closest_3[4], 289.0)

        self.assertEqual(repeats_closest_1, [35, 32, 45, 51, 30, 35, 45])
        self.assertEqual(repeats_closest_2[0], 51)
        self.assertEqual(repeats_closest_2[1], 51)
        self.assertEqual(repeats_closest_2[4], 51)
        self.assertEqual(repeats_closest_2[6], 51)

        self.assertEqual(repeats_closest_3[0], 25)
        self.assertEqual(repeats_closest_3[4], 30)

        self.assertEqual(repeats_1, [35, 33, 51, 46, 30, 35, 42])
        self.assertEqual(repeats_2[0], 61)
        self.assertEqual(repeats_2[1], 54)
        self.assertEqual(repeats_2[4], 88)
        self.assertEqual(repeats_2[6], 55)

        self.assertEqual(repeats_3[0], 21)
        self.assertEqual(repeats_3[4], 30)




    def test_format_columns(self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        triplets_table=find_closest_control_peak_to_sample_peaks(triplets_table,controls)
        triplets_table=get_number_of_triplet_repeats(triplets_table)
        triplets_table=format_columns(triplets_table, controls, "tester1", "gene1")

        size1=list(triplets_table["Size 1"])
        size2=list(triplets_table["Size 2"])
        size3=list(triplets_table["Size 3"])
        repeats_1=list(triplets_table["Repeats_1"])
        repeats_2=list(triplets_table["Repeats_2"])
        repeats_3=list(triplets_table["Repeats_3"])


        self.assertEqual(size1, [306, 300, 342, 362, 288, 306,315])
        self.assertEqual(size2[0], 406.0)
        self.assertEqual(size2[1], 385.0)
        self.assertEqual(size2[4], 486.0)
        self.assertEqual(size2[6], 387.0)
        self.assertEqual(size3[0], 265.0)
        self.assertEqual(size3[4], 289.0)


        self.assertEqual(repeats_1, [35, 33, 51, 46, 30, 35, 42])
        self.assertEqual(repeats_2[0], 61)
        self.assertEqual(repeats_2[1], 54)
        self.assertEqual(repeats_2[4], 88)
        self.assertEqual(repeats_2[6], 55)

        self.assertEqual(repeats_3[0], 21)
        self.assertEqual(repeats_3[4], 30)


if __name__=='__main__':
    unittest.main()
