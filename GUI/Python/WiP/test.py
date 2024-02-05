from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols


molecules = ['COCC','CCCC','CCCF']

ms = [Chem.MolFromSmiles(x) for x in molecules]
print(ms)

fps = [FingerprintMols.FingerprintMol(x) for x in ms]
print(fps)

# the list for the dataframe
qu, ta, sim = [], [], []

# compare all fp pairwise without duplicates
for n in range(len(fps)-1): # -1 so the last fp will not be used
    s = DataStructs.BulkTanimotoSimilarity(fps[n], fps[n+1:]) # +1 compare with the next to the last fp
    print(molecules[n], molecules[n+1:]) # witch mol is compared with what group
    # collect the SMILES and values
    for m in range(len(s)):
        qu.append(molecules[n])
        ta.append(molecules[n+1:][m])
        sim.append(s[m])
