import pandas
import numpy




def get_triplets_table(gene, worksheet):

    '''

    Open the file outputted from genemapper.
    Remove the Normal, Control and NTC.
    Round all the peak columns to the nearest integer.

    '''


    #try opening the triplets file otherwise output that the file could not be found.

    try:
        triplets=pandas.read_csv(worksheet+"_"+gene+".txt", sep="\t")
    except:
        file=open(worksheet+"_"+gene+"_triplets_output.txt", 'w')
        file.write("Genemapper file could not be found- check file name")
        file.close()    



    #extract the peak sizes columns from the table

    triplets_table=triplets.iloc[:,[0,2,6,10]]


    #split the first column to extract the sample id

    sample=triplets_table["Sample File"].str.split("_", n=2, expand=True)
    sample2=list(sample[1])
    triplets_table["Sample"]=sample2
    triplets_table=triplets_table.iloc[:,[0,4,1,2,3]]


    #Remove the Normal, Control and NTC rows

    triplets_table=triplets_table[triplets_table['Sample']!="Normal"] 
    triplets_table=triplets_table[triplets_table['Sample']!="Control"] 
    triplets_table=triplets_table[triplets_table['Sample']!="NTC"] 


    #Round 1st peak column to the nearest integer

    triplets_table['Size 1']=triplets_table['Size 1'].apply(lambda x: round(x))



    #Round 2nd peak column to the nearest integer

    triplets_table_num_rows=triplets_table.shape[0]
    a=0
    while (a<triplets_table_num_rows):
        number=triplets_table.iloc[a,2]
        if numpy.isnan(number):
            number=None
        else:
           number=round(number)
        triplets_table.iloc[a,2]=number
        a=a+1



    #Round 3rd peak column to the nearest integer

    triplets_table_num_rows=triplets_table.shape[0]
    a=0
    while (a<triplets_table_num_rows):
        number=triplets_table.iloc[a,3]
        if numpy.isnan(number):
            number=None
        else:
           number=round(number)
        triplets_table.iloc[a,3]=number
        a=a+1

    return(triplets,triplets_table)






def match_control_samples_with_references(triplets,gene):


    ''' 

    Extract the rows of Normals and Controls from the file outputted from genemapper.
    Match the samples in this table to the ones in the reference controls excel spreadsheet.
    Determine if the values of the controls output from genemapper are within +/- 1 of the values of the controls in the controls table.

    '''


    #Extract only the rows of Normal/control samples from the original triplets table

    sample=triplets["Sample File"].str.split("_", n=3, expand=True)
    sample2=list(sample[1])
    triplets["Sample"]=sample2
    sample3=list(sample[2])
    triplets["Sample2"]=sample3
    controls_1=triplets[triplets['Sample']=="Normal"]
    controls_2=triplets[triplets['Sample']=="Control"]


    #Read in file of reference controls

    controls=pandas.concat([controls_1,controls_2])
    controls=controls.iloc[:,[43,2,6]]

    triplet_control_file=pandas.read_excel("Triplet_controls.xlsx",gene)
    triplet_control_file=pandas.DataFrame(triplet_control_file)
    

    #split the peaks and triplets columns

    peaks=triplet_control_file["Exp_peaks"].str.split("/", n=2, expand=True)
    triplets=triplet_control_file["Exp_repeats"].str.split("/", n=2, expand=True)
    peaks_1=list(peaks[0])
    peaks_2=list(peaks[1])
    triplets_1=list(triplets[0])
    triplets_2=list(triplets[1])

    triplet_control_file["peaks_1"]=peaks_1
    triplet_control_file["peaks_2"]=peaks_2
    triplet_control_file["triplets_1"]=triplets_1
    triplet_control_file["triplets_2"]=triplets_2

    peak_list_1=[]
    peak_list_2=[]
    triplets_list_1=[]
    triplets_list_2=[]


    #only keep the rows of the reference control table that match sample ids of the controls used

    num_rows_controls=controls.shape[0]
    num_rows_triplet_controls=triplet_control_file.shape[0]
    a=0
    while (a<num_rows_controls):
        b=0
        while (b<num_rows_triplet_controls):
            if (controls.iloc[a,0]==triplet_control_file.iloc[b,1]):
                peak_list_1.append(triplet_control_file.iloc[b,7])
                peak_list_2.append(triplet_control_file.iloc[b,8])
                triplets_list_1.append(triplet_control_file.iloc[b,9])
                triplets_list_2.append(triplet_control_file.iloc[b,10])
            b=b+1
        a=a+1
    

    controls["peak_1"]=peak_list_1
    controls["peak_2"]=peak_list_2
   
    controls["triplets_1"]=triplets_list_1
    controls["triplets_2"]=triplets_list_2


    #find out if the control value is within +/- 1 of the reference control

    a=0
    value=0
    continue_list=[]
    num_rows_controls=controls.shape[0]
    while(a<num_rows_controls):
        controls.iloc[a,3]=int(controls.iloc[a,3])
        controls.iloc[a,4]=int(controls.iloc[a,4])
        if(((controls.iloc[a,1]-controls.iloc[a,3])<1 and (controls.iloc[a,1]-controls.iloc[a,3])>0) or ((controls.iloc[a,1]-controls.iloc[a,3])>-1 and (controls.iloc[a,1]-controls.iloc[a,3]<0))):
            continue_list.append("yes")
        else:
            continue_list.append("no")
        a=a+1
     
    controls["valid"]=continue_list     
    
    continue_program="yes"

    filtered_df=controls[controls['valid']=="yes"]
    num_rows_filtered=filtered_df.shape[0]
    if (num_rows_filtered!=num_rows_controls):
        continue_program="no"
    controls=controls.iloc[:,[0,1,2,5,6]]
    return (controls,continue_program)





def find_closest_control_peak_to_sample_peaks(triplets_table,controls):


    '''

    Create a list of the peak values in the control samples.
    Match the peak values of each of the samples to the values in this list to find the closest.
    Add columns for the number of repeats that correspond to the closest peak values.

    '''



    #round the values of the peak columns

    controls['Size 1']=controls['Size 1'].apply(lambda x: round(x))
    controls['Size 2']=controls['Size 2'].apply(lambda x: round(x))


    #make a list of the control sample values to compare the peak sizes of the samples to
    
    list1=list(controls["Size 1"])
    list2=list(controls["Size 2"])
    list3=list(set(list1+list2))

    triplets_table_num_rows=triplets_table.shape[0]
    numbers=[]
    a=0
    while (a<triplets_table_num_rows):
        number=triplets_table.iloc[a,2]
        numbers.append(min(list3, key=lambda x:abs(x-number)))
    
        a=a+1

    triplets_table["closest_1"]=numbers

  
    a=0
    numbers=[]
    num_rows_triplets_table=triplets_table.shape[0]
    while (a<num_rows_triplets_table):
        number=triplets_table.iloc[a,3]
        if (numpy.isnan(number)):
            numbers.append("NaN")
        else:   
            numbers.append(min(list3, key=lambda x:abs(x-number)))
   
        a=a+1
    
    triplets_table["closest_2"]=numbers


    a=0
    numbers=[]
    num_rows_triplets_table=triplets_table.shape[0]
    while (a<num_rows_triplets_table):
        number=triplets_table.iloc[a,4]
        if (numpy.isnan(number)):
            numbers.append("NaN")
        else:   
            numbers.append(min(list3, key=lambda x:abs(x-number)))
   
        a=a+1
    
    triplets_table["closest_3"]=numbers




    #Find the closest control peak size to each of the 1st sample peak sizes

    triplets_table_num_rows=triplets_table.shape[0]
    a=0
    repeats_closest1=[]
    num_rows_controls=controls.shape[0]
    while (a<triplets_table_num_rows):
        b=0
        while (b<num_rows_controls):
            if (triplets_table.iloc[a,5]==controls.iloc[b,1]):
                repeats_closest1.append(controls.iloc[b,3])
            elif (triplets_table.iloc[a,5]==controls.iloc[b,2]):
                repeats_closest1.append(controls.iloc[b,4])
            b=b+1
        a=a+1
        
    triplets_table["repeats_closest_1"]=repeats_closest1



    #Find the closest control peak size to each of the 2nd sample peak sizes

    a=0

    repeats_closest2=[]
    while (a<triplets_table_num_rows):
        b=0
        c=0
        while (b<num_rows_controls):
            if (triplets_table.iloc[a,6]==controls.iloc[b,1]):
                repeats_closest2.append(controls.iloc[b,3])
                c=1
            elif (triplets_table.iloc[a,6]==controls.iloc[b,2]):
                repeats_closest2.append(controls.iloc[b,4])
                c=1
            b=b+1
        if (c==0):
            repeats_closest2.append("NaN")

        a=a+1
    
    
    triplets_table["repeats_closest_2"]=repeats_closest2



    #Find the closest control peak size to each of the 3rd sample peak sizes

    a=0

    repeats_closest3=[]
    while (a<triplets_table_num_rows):
        b=0
        c=0
        while (b<num_rows_controls):
            if (triplets_table.iloc[a,7]==controls.iloc[b,1]):
                repeats_closest3.append(controls.iloc[b,3])
                c=1
            elif (triplets_table.iloc[a,7]==controls.iloc[b,2]):
                repeats_closest3.append(controls.iloc[b,4])
                c=1
            b=b+1
        if (c==0):
            repeats_closest3.append("NaN")

        a=a+1
    
    
    triplets_table["repeats_closest_3"]=repeats_closest3

    return(triplets_table)





def get_number_of_triplet_repeats(triplets_table):


    '''

    Find the difference between the sample peak size and the peak size of the closest control.
    Divide this value by 3 and add it to the number of repeats in the control, to find the number of repeats the sample peak correlates to.
    Repeat this for all three peaks for all samples.

    '''    

 
    #Find the difference between the sample peak size and nearest control peak size and divide this by 3. Add this value to the number of repeats that correlate to the closest control peak size value

    a=0

    difference1=[]
    triplets_table_num_rows=triplets_table.shape[0]
    while (a<triplets_table_num_rows):
        difference=triplets_table.iloc[a,2]-triplets_table.iloc[a,5]
        if (difference==0):
            triplets_table.iloc[a,8]=int(triplets_table.iloc[a,8])
            difference=triplets_table.iloc[a,8]+difference
            difference=difference.round()
            difference1.append(difference)
        else:
            difference=difference/3
            triplets_table.iloc[a,8]=int(triplets_table.iloc[a,8])
            difference=triplets_table.iloc[a,8]+difference
            difference=difference.round()
            difference=int(difference)
            difference1.append(difference)
        
        a=a+1
    
    
    triplets_table["Repeats_1"]=difference1


    #Do the same for the 2nd sample peak size column

    a=0

    difference2=[]
    while (a<triplets_table_num_rows):
        if (numpy.isnan(triplets_table.iloc[a,3])):
            difference2.append("NaN")
        else:
            triplets_table.iloc[a,3]=int(triplets_table.iloc[a,3])
            triplets_table.iloc[a,6]=int(triplets_table.iloc[a,6])
            difference=triplets_table.iloc[a,3]-triplets_table.iloc[a,6]
            if (difference==0):
                triplets_table.iloc[a,9]=int(triplets_table.iloc[a,9])
                difference=triplets_table.iloc[a,9]
                difference2.append(difference)
            else:
                difference=difference/3
                triplets_table.iloc[a,9]=int(triplets_table.iloc[a,9])
                difference=triplets_table.iloc[a,9]+difference
                difference=difference.round()
                difference=int(difference)
                difference2.append(difference)
        
        a=a+1
    
 
    triplets_table["Repeats_2"]=difference2



   #Do the same for the 3rd sample peak size column

    a=0

    difference3=[]
    triplets_table_num_rows=triplets_table.shape[0]
    while (a<triplets_table_num_rows):
        if (numpy.isnan(triplets_table.iloc[a,4])):
            difference3.append("NaN")
        else:
            triplets_table.iloc[a,4]=int(triplets_table.iloc[a,4])
            triplets_table.iloc[a,7]=int(triplets_table.iloc[a,7])
            difference=triplets_table.iloc[a,4]-triplets_table.iloc[a,7]
            if (difference==0):
                triplets_table.iloc[a,10]=int(triplets_table.iloc[a,10])
                difference=triplets_table.iloc[a,10]
                difference3.append(difference)
            else:
                difference=difference/3
                triplets_table.iloc[a,10]=int(triplets_table.iloc[a,10])
                difference=triplets_table.iloc[a,10]+difference
                difference=difference.round()
                difference=int(difference)
                difference3.append(difference)
        
        a=a+1
    
 
    triplets_table["Repeats_3"]=difference3

    return(triplets_table)





def format_columns(triplets_table, controls, worksheet, gene):


    '''

    Extract the sample, peak sizes and repeat columns that will be in the final output.
    Replace empty values with NaN.
    Output table as a file in the current folder.

    '''


    #Extract the sample, peak sizes and repeats columns
    triplets_table=triplets_table.iloc[:,[0,2,3,4,11,12,13]]


    a=0
    numbers=[]
    triplets_table_num_rows=triplets_table.shape[0]
    while (a<triplets_table_num_rows):
        number=triplets_table.iloc[a,2]
        if  (numpy.isnan(number)):
            number="NaN"
        else:
            number=int(number)
        numbers.append(number)
        a=a+1
    triplets_table["Size 2"]=numbers


    a=0
    numbers=[]
    triplets_table_num_rows=triplets_table.shape[0]
    while (a<triplets_table_num_rows):
        number=triplets_table.iloc[a,3]
        if  (numpy.isnan(number)):
            number="NaN"
        else:
            number=int(number)
        numbers.append(number)
        a=a+1
    triplets_table["Size 3"]=numbers

    triplets_table.to_csv(worksheet+"_"+gene+"_triplets_output.txt", index=None, sep='\t')

    return(triplets_table)




if __name__ == "__main__":


    #The user can enter the gene name and the worksheet number 
    
    gene=input('Enter gene:')
    worksheet=input('Enter worksheet:')


    #Check the gene name entered is one of the the ones in the triplet_controls file


    if ((gene=="FRAX") or (gene=="FA") or (gene=="C9ORF72") or (gene=="HD") or (gene=="MDMYo(DM1)") or (gene=="SCA1") or (gene=="SCA2") or (gene=="SCA3") or (gene=="SCA6")):


        triplets,triplets_table=get_triplets_table(gene, worksheet)

        controls,continue_program=match_control_samples_with_references(triplets,gene)


 
        #Continue program only if control values in genemapper output are within +/- 1 of the values in the controls file


        if (continue_program=="yes"):

            triplets_table_2=find_closest_control_peak_to_sample_peaks(triplets_table,controls)

            triplets_table_3=get_number_of_triplet_repeats(triplets_table_2)

            format_columns(triplets_table_3, controls, worksheet, gene)

        else:
            file=open(worksheet+"_"+gene+"_triplets_output.txt", 'w')
            file.write("The controls are not within +/- 1 of the reference controls")
            file.close()

    else:
        file=open(worksheet+"_"+gene+"_triplets_output.txt",'w')
       	file.write("Gene entered incorrectly")
       	file.close()
