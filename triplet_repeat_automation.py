
'''
triplet_repeat_automation.py 

Calculates the number of triplet repeats for sample peak sizes outputted from Genemapper.
Author: Laura McCluskey
Version 1.0.0
'''

import pandas
import numpy
import xlwings
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

def get_triplets_table(gene, worksheet):

    '''
    Open the file outputted from genemapper.
    Remove the Normal, Control and NTC.
    Round all the peak columns to the nearest integer.

    '''

    #try opening the triplets file, otherwise output that the file could not be found.
    try:
        triplets=pandas.read_csv(worksheet+"_"+gene+".txt", sep="\t")
    except:
        file=open(worksheet+"_"+gene+"_triplets_output.txt", 'w')
        file.write("Genemapper file could not be found- check file name")
        file.close()    

    #check the extra column needed on the end of the table hasn't been deleted in editing process
    if (len(triplets.columns)==41):
        triplets["Extra_column"]=""
    
    #extract the peak sizes columns from the table
    triplets_table=triplets.filter(items=["Sample File", "Size 1", "Size 2", "Size 3"])

    #split the first column to extract the sample id
    sample=triplets_table["Sample File"].str.split("_", n=2, expand=True)
    sample2=list(sample[1])
    triplets_table["Sample"]=sample2
    triplets_table=triplets_table.filter(items=["Sample File","Sample","Size 1","Size 2","Size 3"])

    #Remove the Normal, Control and NTC rows
    triplets_table=triplets_table[triplets_table['Sample']!="Normal"] 
    triplets_table=triplets_table[triplets_table['Sample']!="Control"] 
    triplets_table=triplets_table[triplets_table['Sample']!="NTC"] 

    #Round peak columns to the nearest integer
    triplets_table['Size 1']=triplets_table['Size 1'].apply(lambda x: None if numpy.isnan(x) else round(x) )
    triplets_table['Size 2']=triplets_table['Size 2'].apply(lambda x: None if numpy.isnan(x) else round(x) )
    triplets_table['Size 3']=triplets_table['Size 3'].apply(lambda x: None if numpy.isnan(x) else round(x) )

    return(triplets,triplets_table)

def match_control_samples_with_references(triplets,gene):

    '''
     Extract the rows of Normals and Controls from the file outputted from genemapper.
     Match the samples in this table to the ones in the reference controls excel spreadsheet.
     Determine if the values of the controls output from genemapper are within +/- 3 of the values of the controls in the controls table.
     '''

    #Add sampleid column to triplets dataframe
    sample=triplets["Sample File"].str.split("_", n=3, expand=True)
    sample2=list(sample[1])
    triplets["Sample"]=sample2
    sample3=list(sample[2])
    triplets["Sample2"]=sample3


    #Extract only the rows of Normal/control samples from the triplets dataframe
    controls_1=triplets[triplets['Sample']=="Normal"]
    controls_2=triplets[triplets['Sample']=="Control"]

    controls=pandas.concat([controls_1, controls_2])
    controls=controls.filter(items=["Sample2", "Size 1", "Size 2"])

    #Read in file of reference controls- xlwings used to allow reading of password protected excel spreadsheet
    if (gene=="FRAX"):
        triplets_excel_input=xlwings.Book("Triplet_controls_FRAX.xlsx")
    else:
        triplets_excel_input=xlwings.Book("Triplet_controls.xlsx")

    triplet_control_file=triplets_excel_input.sheets[gene]
    triplet_control_file=triplet_control_file['A1:G1000000'].options(pandas.DataFrame,index=False, header=True).value


    #split the peaks and triplets columns
    peaks=triplet_control_file["Exp_peaks"].str.split("/", n=2, expand=True)
    triplets=triplet_control_file["Exp_repeats"].str.split("/", n=2, expand=True)
    triplet_control_file["peaks_1"]=list(peaks[0])
    triplet_control_file["peaks_2"]=list(peaks[1])
    triplet_control_file["triplets_1"]=list(triplets[0])
    triplet_control_file["triplets_2"]=list(triplets[1])


    #only keep the rows of the reference control table that match sample ids of the controls used
    new_table=pandas.merge(left=controls, right=triplet_control_file, how='left', left_on='Sample2', right_on='reference_sample')
    controls=new_table.filter(items=["Sample2", "Size 1","Size 2", "peaks_1", "peaks_2", "triplets_1", "triplets_2"])

    #find out if the control value is within +/- 3 of the reference control
    controls['peaks_1']=controls['peaks_1'].apply(lambda x: int(x))
    controls['peaks_2']=controls['peaks_2'].apply(lambda x: int(x))
    controls["difference in peak 1"]= controls["peaks_1"]-controls["Size 1"]
    controls["difference in peak 2"]= controls["peaks_2"]-controls["Size 2"]
    controls_filtered=controls[controls["difference in peak 1"].between(-3,3, inclusive=True)]
    controls_filtered=controls_filtered[controls_filtered["difference in peak 2"].between(-3,3, inclusive=True)]

    #only continue with program if control value is within +/- 3 of the reference control
    if (controls_filtered.shape[0] ==controls.shape[0]):
        continue_program="yes"
    else:
        continue_program="no"

    controls=controls.filter(items=["Sample2", "Size 1","Size 2","triplets_1", "triplets_2"])
    return (controls,continue_program)


def get_closest_value(x, array):

    '''
    Function called from find_closest_control_peak_to_sample_peaks function
    Input: Each value of Size columns in triplets table and an array of all peak sizes in controls table
    Output: size of peak in controls table closest to value of input
    '''
    if (x==numpy.nan):
        value=numpy.nan
    else:
        array_minus_x=abs(array-x)
        array_minus_x=array_minus_x.tolist()
        min_index=array_minus_x.index(min(array_minus_x))
        value=array[min_index]
    return(value)


def find_closest_control_peak_to_sample_peaks(triplets_table,controls):

    '''
    Create a list of the peak values in the control samples.
    Match the peak values of each of the samples to the values in this list to find the closest.
    Add columns for the number of repeats that correspond to the closest peak values.
    '''

    #round the values of the peak columns
    controls['Size 1']=controls['Size 1'].apply(lambda x: round(x))
    controls['Size 2']=controls['Size 2'].apply(lambda x: round(x))

    #make a list of the peak sizes in the reference controls file  
    list1=list(controls["Size 1"])
    list2=list(controls["Size 2"])
    list3=list(set(list1+list2))
    peak_array=numpy.array(list3)

    #match the peak sizes in the triplets table and the reference controls file using get_closest_value_function
    triplets_table["Size 1"]=triplets_table['Size 1'].apply(lambda x: numpy.nan if x==None else x)
    triplets_table["Size 2"]=triplets_table['Size 2'].apply(lambda x: numpy.nan if x==None else x)
    triplets_table["Size 3"]=triplets_table['Size 3'].apply(lambda x: numpy.nan if x==None else x)

    triplets_table["closest_1"]=triplets_table['Size 1'].apply(lambda x: get_closest_value(x,peak_array))
    triplets_table["closest_2"]=triplets_table['Size 2'].apply(lambda x: get_closest_value(x,peak_array))
    triplets_table["closest_3"]=triplets_table['Size 3'].apply(lambda x: get_closest_value(x,peak_array))

    #Make a table with two columns (Size and triplets) from the reference controls table
    controls1=controls.filter(items=["Size 1","triplets_1"])
    controls1.columns=["Size","triplets"]
    controls2=controls.filter(items=["Size 2","triplets_2"])
    controls2.columns=["Size","triplets"]
    controls_altered=pandas.concat([controls1, controls2])

    #Merge the controls table with triplets table, joining on closest_2 column to get the number of triplets each peak size correlates to 
    triplets_table=pandas.merge(left=triplets_table, right=controls_altered, how='left', left_on='closest_1', right_on='Size')
    triplets_table=triplets_table.filter(items=["Sample File","Sample", "Size 1", "Size 2", "Size 3", "closest_1", "closest_2", "closest_3", "triplets"])
    triplets_table.columns=["Sample File","Sample", "Size 1", "Size 2", "Size 3", "closest_1", "closest_2", "closest_3", "repeats_closest_1"]

    #Merge the controls table with triplets table, joining on closest_2 column to get the number of triplets each peak size correlates to 
    triplets_table=pandas.merge(left=triplets_table, right=controls_altered, how='left', left_on='closest_2', right_on='Size')
    triplets_table=triplets_table.filter(items=["Sample File","Sample", "Size 1", "Size 2", "Size 3", "closest_1", "closest_2", "closest_3","repeats_closest_1", "triplets"])
    triplets_table.columns=["Sample File","Sample", "Size 1", "Size 2", "Size 3", "closest_1", "closest_2", "closest_3", "repeats_closest_1", "repeats_closest_2"]

    #Merge the controls table with triplets table, joining on closest_3 column to get the number of triplets each peak size correlates to 
    triplets_table=pandas.merge(left=triplets_table, right=controls_altered, how='left', left_on='closest_3', right_on='Size')
    triplets_table=triplets_table.filter(items=["Sample File","Sample", "Size 1", "Size 2", "Size 3", "closest_1", "closest_2", "closest_3","repeats_closest_1", "repeats_closest_2","triplets"])
    triplets_table.columns=["Sample File","Sample", "Size 1", "Size 2", "Size 3", "closest_1", "closest_2", "closest_3", "repeats_closest_1", "repeats_closest_2", "repeats_closest_3"]

    return(triplets_table)


def get_number_of_triplet_repeats(triplets_table):

    '''
    Find the difference between the sample peak size and the peak size of the closest control.
    Divide this value by 3 to get the difference in the number of triplet repeats.
    Add his difference to the number of repeats in the control, to find the number of repeats the sample peak correlates to.
    Repeat this for all three peaks for all samples.
    '''

    triplets_table['Size 1']=triplets_table['Size 1'].apply(lambda x: int(x)) 
    triplets_table['closest_1']=triplets_table['closest_1'].apply(lambda x: int(x))
    triplets_table['repeats_closest_1']=triplets_table['repeats_closest_1'].apply(lambda x: int(x))
    triplets_table['difference']=(triplets_table["Size 1"]-triplets_table["closest_1"])/3
    triplets_table["Repeats_1"]=triplets_table["repeats_closest_1"]+triplets_table["difference"]
    triplets_table['Repeats_1']=triplets_table['Repeats_1'].apply(lambda x: round(x))
    triplets_table=triplets_table.filter(items=["Sample File","Sample", "Size 1", "Size 2", "Size 3", "closest_1", "closest_2", "closest_3","repeats_closest_1", "repeats_closest_2", "repeats_closest_3", "Repeats_1"])

    triplets_table['Size 2']=triplets_table['Size 2'].apply(lambda x: numpy.nan if numpy.isnan(x) else int(x)) 
    triplets_table["closest_2"]=triplets_table['closest_2'].apply(lambda x:  numpy.nan if x=="NaN" else int(x))
    triplets_table["repeats_closest_2"]=triplets_table['repeats_closest_2'].apply(lambda x:  float(x))
    triplets_table['difference']=(triplets_table["Size 2"]-triplets_table["closest_2"])/3
    triplets_table["Repeats_2"]=triplets_table["repeats_closest_2"]+triplets_table["difference"]
    triplets_table['Repeats_2']=triplets_table['Repeats_2'].apply(lambda x: "NaN" if numpy.isnan(x) else round(x))
    triplets_table=triplets_table.filter(items=["Sample File","Sample", "Size 1", "Size 2", "Size 3", "closest_1", "closest_2", "closest_3","repeats_closest_1", "repeats_closest_2", "repeats_closest_3", "Repeats_1", "Repeats_2"])

    triplets_table['Size 3']=triplets_table['Size 3'].apply(lambda x: numpy.nan if numpy.nan else int(x)) 
    triplets_table["closest_3"]=triplets_table['closest_3'].apply(lambda x:  numpy.nan if x=="NaN" else int(x))
    triplets_table["repeats_closest_3"]=triplets_table['repeats_closest_3'].apply(lambda x:  float(x))
    triplets_table['difference']=(triplets_table["Size 3"]-triplets_table["closest_3"])/3
    triplets_table["Repeats_3"]=triplets_table["repeats_closest_3"]+triplets_table["difference"]
    triplets_table['Repeats_3']=triplets_table['Repeats_3'].apply(lambda x: "NaN" if numpy.isnan(x) else round(x))
    triplets_table=triplets_table.filter(items=["Sample File","Sample", "Size 1", "Size 2", "Size 3", "closest_1", "closest_2", "closest_3","repeats_closest_1", "repeats_closest_2", "repeats_closest_3", "Repeats_1", "Repeats_2", "Repeats_3"])

    return(triplets_table)


def format_columns(triplets_table, controls, worksheet, gene):


    #Extract the sample, peak sizes and repeats columns and output table to a text file
    triplets_table=triplets_table.filter(items=["Sample File", "Size 1", "Size 2", "Size 3", "Repeats_1", "Repeats_2", "Repeats_3"])
    triplets_table['Size 1']=triplets_table['Size 1'].apply(lambda x: "NaN" if numpy.isnan(x) else int(x))
    triplets_table['Size 2']=triplets_table['Size 2'].apply(lambda x: "NaN" if numpy.isnan(x) else int(x))
    triplets_table['Size 3']=triplets_table['Size 3'].apply(lambda x: "NaN" if numpy.isnan(x) else int(x))
    triplets_table.to_csv(worksheet+"_"+gene+"_triplets_output.txt", index=None, sep='\t')

    #output the results to an excel spreadsheet with checking boxes
    wb=Workbook()
    ws1=wb.create_sheet("Triplet_results")
    for row in dataframe_to_rows(triplets_table):
        ws1.append(row)
 
    ws1["K2"]="Worksheet:"
    ws1["L2"]=worksheet
    ws1["K5"]="First checker:"
    ws1["K6"]="Date:"
    ws1["K7"]="Second checker:"
    ws1["K8"]="Date:"    

    ws1.column_dimensions['B'].width=60
    ws1.column_dimensions['F'].width=10
    ws1.column_dimensions['G'].width=10
    ws1.column_dimensions['H'].width=10
    ws1.column_dimensions['K'].width=20

    wb.save(worksheet+"_"+gene+"_triplets_output_excel.xlsx")

    return(triplets_table)


if __name__ == "__main__":

    #The user can enter the gene name and the worksheet number 
    
    gene=input('Enter gene:')
    worksheet=input('Enter worksheet:')

    #Check the gene name entered is one of the the ones in the triplet_controls file

    if ((gene=="FRAX") or (gene=="FA") or (gene=="C9ORF72") or (gene=="HD") or (gene=="MDMYo(DM1)") or (gene=="SCA1") or (gene=="SCA2") or (gene=="SCA3") or (gene=="SCA6") or (gene=="gene1")):

        triplets,triplets_table=get_triplets_table(gene, worksheet)

        controls,continue_program=match_control_samples_with_references(triplets,gene)

        #Continue program only if control values in genemapper output are within +/- 1 of the values in the controls file
        if (continue_program=="yes"):

            triplets_table_2=find_closest_control_peak_to_sample_peaks(triplets_table,controls)

            triplets_table_3=get_number_of_triplet_repeats(triplets_table_2)

            format_columns(triplets_table_3, controls, worksheet, gene)

        else:
            file=open(worksheet+"_"+gene+"_triplets_output.txt", 'w')
            file.write("The controls are not within +/- 3 of the reference controls")
            file.close()

    else:
        file=open(worksheet+"_"+gene+"_triplets_output.txt",'w')
       	file.write("Gene entered incorrectly")
       	file.close()
