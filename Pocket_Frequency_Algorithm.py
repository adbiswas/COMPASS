import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os
import pandas as pd
import numpy as np
import csv
import gc

###
entries = 'INPUT_YOUR_PATH_OF_THE_DIRECTORY_WHERE_YOU HAVE_THE_POCKETS'

dfinal = pd.DataFrame(columns = ['ResID', 'Target'])
data =[]

#### READs all the FILES ######
for root, directories, files in os.walk(entries, topdown=True):
    for name in files:
        df = pd.read_excel(os.path.join(root,name))
        basename_without_ext = os.path.splitext(os.path.basename(name))[0]
        df.loc[:,"Pocket"] = basename_without_ext

        basename_without = os.path.splitext(os.path.basename(root))[0]
        df.loc[:,"Source"] = basename_without
        
        ############ OPEN THE POCKET FILES #################
        df.columns = ['Number','Atom', 'Element', 'Residue', 'ResNum', 'Chain', 'X', 'Y', 'Z', 'Charge', 'Potential', 'Constraint', "Score", "Consensus", "Pocket", "Source"]
        
        ############ SELECT THE  RESIDUES BASED ON THE C-alpha ATOMS #################
        d1=df.loc[df.Atom == 'CA']
        index = d1.index
        number_of_rows = len(index)
        
        ########### JOIN INFORMATION OF RESIDUE, CHAIN and POCKET NUMBER ##################
        d2 = d1[['Residue', 'ResNum', 'Chain', 'Pocket', 'Source', "Score", "Consensus"]]
        d3 = d2['Residue'].str.cat(d2['ResNum'].values.astype(str)).str.cat(d2['Chain'])
        d4 = d2['Pocket'].str.cat(d2['Chain']).str.cat(d2['Source'], sep ="-")
        d3a = d3.to_frame()
        d4a = d4.to_frame()
        
        d3a.columns =['ResID']
        d4a.columns = ['Target']
        
        d5 = pd.concat([d3a, d4a, d2["Score"], d2["Consensus"]], axis=1, join='inner')
        dfinal =dfinal.append(d5) #Create a list with Residue number and Chain Name
gc.collect()
############### Select residues with equal residues number and equal chain name ############## 
for i in range(len(dfinal)):
    k=0
    for j in range(len(dfinal)):
        if dfinal.iloc[i, 0] == dfinal.iloc[j, 0]:
            if dfinal.iloc[i, 1] != dfinal.iloc[j, 1]:
                k = k+1
                d6 = [dfinal.iloc[i, 0], dfinal.iloc[i, 1], dfinal.iloc[j, 0], dfinal.iloc[j, 1], k, dfinal.iloc[i, 2], dfinal.iloc[i, 3]]
                data.append(d6) ## Save data: similar resid, number of times appeared, docking score and consensus           
gc.collect()
############# Save the above file in the file result.csv ##############################

file = open('result.csv', 'w+', newline ='')
with file:    
    write = csv.writer(file)
    write.writerows(data)

########################################################################################################
####################################### POCKET SEARCH ##################################################

poc1df = pd.read_csv('result.csv', sep='Delim_first|Delim_second|[,]', engine='python', header=None)

dat2 = []

###################### Remove repetition of pocket IDs ##########################################
for i in range(len(poc1df)):
    for j in range(len(poc1df)):
        if poc1df.iloc[i, 1] != poc1df.iloc[j, 1]:
            dat1 = [poc1df.iloc[i, 1], poc1df.iloc[i, 3], poc1df.iloc[i, 5], poc1df.iloc[i, 6]]
            #dat2.append(dat1)
            
    dat2.append(dat1)
    #print(sys.getsizeof(dat2))
    #print('First DAT', dat2)
    dat3 = pd.DataFrame(dat2)
    dat4 = dat3.drop_duplicates()
    dat2 = dat4.values.tolist()
    #print('Second DAT', dat2)
    #print(sys.getsizeof(dat2))
    gc.collect()

'''
#print(type(dat2))
dat3 = pd.DataFrame(dat2)
dat4 = dat3.drop_duplicates()
dat2 = dat4.values.tolist()
#print(type(dat3))
'''
################################################################################################

################# Save the data in the file Poc-match-list1 ###################################

file = open('Poc-match-list1.csv', 'w+', newline ='')
with file:    
    write = csv.writer(file)
    write.writerows(dat2)
#################################################################################################
#gc.collect()
################# Change the line and insert * to differentiate #################################
colnames=['A', 'B', 'C', 'D'] 
poc2df = pd.read_csv('Poc-match-list1.csv', sep='Delim_first|Delim_second|[,]', names=colnames, engine='python')

poc3df = poc2df['A'].str.cat(poc2df['B'], sep ="*").str.cat(poc2df['C'].values.astype(str), sep ="*").str.cat(poc2df['D'].values.astype(str), sep ="*")
poc4df = poc3df.to_frame()
poc5df = poc4df.drop_duplicates()
poc5df.to_csv('Poc-match-list2.csv', index=False, header=False)

################# Load the file into differentiating the Pocket and PDBID #########################
################# Pocket Information consist of Pocket number and Chain Name ##################################
poc6df = pd.read_csv('Poc-match-list2.csv', sep='Delim_first|Delim_second|[*]|Delim_third|[-]', engine='python', names=['A', 'B', 'C', 'D', 'E', 'F'])
###############################################################################################################

############################## Save the above data in Poc-match-list2.csv ##############################################
daat1 = [] #Blank dataset
daat2 = []

##################### Sort the pocket which appear in the same protein ################################################
for i in range(len(poc6df)):
    if poc6df.iloc[i, 1] == poc6df.iloc[i, 3]:
        daat = [poc6df.iloc[i, 0], poc6df.iloc[i, 1], poc6df.iloc[i, 2], poc6df.iloc[i, 3], poc6df.iloc[i, 4], poc6df.iloc[i, 5]]
        daat1.append(daat)
    
##################### Sort the pocket which appear in different protein ################################################
    if poc6df.iloc[i, 1] != poc6df.iloc[i, 3]:
        daat = [poc6df.iloc[i, 0], poc6df.iloc[i, 1], poc6df.iloc[i, 2], poc6df.iloc[i, 3], poc6df.iloc[i, 4], poc6df.iloc[i, 5]]
        daat2.append(daat) 

        
#################### Save the pocket informations ######################################################################

file = open('Similar-Pockets-Same-Protein.csv', 'w+', newline ='')
with file:
    write = csv.writer(file)
    write.writerows(daat1)
    
    
file = open('Similar-Pockets-Different-Proteins.csv', 'w+', newline ='')
with file:
    write = csv.writer(file)
    write.writerows(daat2)

##################################### SIMILAR CAVITIES WITHIN SAME PROTEINS #######################
##################################### Drop Duplicates #############################################
SPSP1 = pd.read_csv('Similar-Pockets-Same-Protein.csv', sep='Delim_first|Delim_second|[,]', engine='python', names=['A', 'B', 'C', 'D', 'E', 'F'])
SPSP1.A, SPSP1.C = np.minimum(SPSP1.A, SPSP1.C), np.maximum(SPSP1.A, SPSP1.C)
SPSP1DF = SPSP1.drop_duplicates()

SPSP1DF.to_csv('Similar-Pockets-Same-Protein-Final.csv', index=False, header=['Pocket', 'Source', 'Similar-Pocket', 'Source', 'Score', 'Consensus'])

########################################################################################################

##################################### SIMILAR CAVITIES WITHIN DIFFERENT PROTEINS #######################

SPDP2 = pd.read_csv('Similar-Pockets-Different-Proteins.csv', sep='Delim_first|Delim_second|[,]', engine='python', names=['A', 'B', 'C', 'D', 'E', 'F'])

SPDP2 = SPDP2.astype(str)

SPDP2[SPDP2.apply(sorted, axis=1).map(str).duplicated(keep='last')]

SPDP3 = SPDP2[SPDP2[['A', 'C']].apply(sorted, axis=1).map(str).duplicated()]

#print(SPDP3)
SPDP3.to_csv('Similar-Pockets-Different-Proteins-Repeat.csv', index=False, header=['Pocket', 'Source', 'Similar-Pocket', 'Source', 'Score', 'Consensus'])

############################################################################################################
 
################################### REMOVE DUPLICATES FURTHER ##############################################

SPDP4 = pd.read_csv('Similar-Pockets-Different-Proteins-Repeat.csv', sep='Delim_first|Delim_second|[,]', engine='python', names=['A', 'B', 'C', 'D', 'E', 'F'])

SPDP4a = SPDP4[['A', 'B', 'C', 'D']]
SPDP5 = SPDP4a.loc[
    pd.DataFrame(
        np.sort(SPDP4a.values, axis=1), columns=SPDP4a.columns
    ).drop_duplicates().index
]
#print(SPDP5)

testa = SPDP4['A'].str.cat(SPDP4['B'], sep ="*").str.cat(SPDP4['C'], sep ="*").str.cat(SPDP4['D'], sep ="*")
testb = SPDP5['A'].str.cat(SPDP5['B'], sep ="*").str.cat(SPDP5['C'], sep ="*").str.cat(SPDP5['D'], sep ="*")

#print(testa)
#print(testb)

NODUP1=[]
for i in range (len(testb)):
    for j in range(len(testa)):
        if testb.iloc[i] == testa.iloc[j]:
            #print(i,j,testb.iloc[i],testa.iloc[j],'I am there')
            #print(SPDP4.iloc[j, 0], SPDP4.iloc[j, 1], SPDP4.iloc[j, 2], SPDP4.iloc[j, 3],  SPDP4.iloc[j, 4],  SPDP4.iloc[j, 5])
            NODUP2 = [SPDP4.iloc[j, 0], SPDP4.iloc[j, 1], SPDP4.iloc[j, 2], SPDP4.iloc[j, 3],  SPDP4.iloc[j, 4],  SPDP4.iloc[j, 5]]
            NODUP1.append(NODUP2)

file = open('Similar-Pockets-Different-Proteins-Final.csv', 'w+', newline ='')
with file:    
    write = csv.writer(file)
    write.writerows(NODUP1)

#SPDP5.to_csv('Similar-Pockets-Different-Proteins-Final.csv', index=False, header=None)

##########################################################################################################

######################## SCORING THE SIMILAR POCKETS IN DIFFERENT PROTEINS ###############################

SPDP6 = pd.read_csv('Similar-Pockets-Different-Proteins-Final.csv', sep='Delim_first|Delim_second|[,]', engine='python', names=['A', 'B', 'C', 'D', 'E', 'F'])

SimilarPocket = SPDP6['A'].str.cat(SPDP6['B'].values.astype(str), sep ="-")

SPTDF = SimilarPocket.to_frame()
#print(len(SimilarPocket))

######## SPTDF = Similar PockeT in Data Frame #############################################
######## SPTDF has information on Pocket number, Chain name and PDB ID ####################

######## Scoring has been done based on how many times did the same pocket with the same chain name and same pdb id #######
######## appears on different pockets, redundancy removed #################################################################

pocdata=[]
for i in range (len(SPTDF)):
    k = 0
    for j in range(len(SPTDF)-1):
        if SPTDF.iloc[i, 0] == SPTDF.iloc[j+1, 0]:
            k=k+1
    datapoc = [SPTDF.iloc[i, 0], k, SPDP6.iloc[i, 4],  SPDP6.iloc[i, 5]]
    pocdata.append(datapoc)

###########################################################################################################################

file = open('Pocket-Score.csv', 'w+', newline ='')
with file:    
    write = csv.writer(file)
    write.writerows(pocdata)
    
    
file.close()


FinalData = pd.read_csv('Pocket-Score.csv', sep='Delim_first|Delim_second|[,]', engine='python', header=None)

PocketF = []

for i in range(len(FinalData)-1):
    if FinalData.iloc[i, 0] != FinalData.iloc[i+1, 0]:
        FD = [FinalData.iloc[i, 0], FinalData.iloc[i, 1], FinalData.iloc[i, 2], FinalData.iloc[i, 3]]
        PocketF.append(FD)


file = open('SPDP-Score.csv', 'w+', newline ='')
with file:    
    write = csv.writer(file)
    write.writerows(PocketF)
    
    
file.close()

##########################################################################################################

print('Done')
