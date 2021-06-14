from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

molecules = open('glucose_degradation_output.csv','r')
lines = molecules.readlines()
counter = 0
with open('Glucose_Desc.csv', 'w') as the_file:
	the_file.write("Generation,Id,NumStereoIsomers"+'\n') 	
	for line in lines:
		counter +=1
		line=line.rstrip('\n')
		line=line.split('\t')
		m = Chem.MolFromSmiles(line[1])
		isomers = tuple(EnumerateStereoisomers(m))
		numste = str(len(isomers))
		the_file.write(line[0]+","+line[1]+","+numste+'\n') 
