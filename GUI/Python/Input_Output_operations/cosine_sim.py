import os, time, csv, shutil
import numpy as np
from natsort import natsorted
from datetime import datetime
from Python.atom_mass import atom_masses
from PyQt6.QtWidgets import QApplication

def cosine_sim(directory_or_csv, similarity_threshold, write_to_file):
    
    '''
    Calculate cosine similarity between molecular structures in a directory or specified in a CSV file.

    Parameters:
    - directory_or_csv: Path to a directory containing .gjf, .inp. and/or .xyz files or a CSV file with a list of file paths.
    - similarity_threshold: Threshold for considering molecules as similar.
    - write_to_file: Boolean flag to indicate whether to write results to a file.
    '''

    def read_file(file_path):
        '''Read the xyz data of a .gjf, .xyz, or .inp file and write the data to a list.'''
        
        #empty list to write coordinate data to
        coordinates = []

        with open(file_path, 'r') as file:
            lines = file.readlines()

        for line in lines:
            line = line.strip()
            split_line = line.split() #split by whitespace
            
            #skip any empty lines
            if split_line:           
                #the only entry with 4 splits, where instance 1 is an atom symbol (sometimes an atom number!), and a period in instances 1-3 will be the xyz coordiante lines. Write these to the geom list
                if len(split_line) == 4 and all('.' in x for x in split_line[1:4]) and (split_line[0].isalpha() or split_line[0].isdigit()):
                    coordinates.append(line.split())
        
        #If no coordinates are found, return None
        if len(coordinates) == 0:
            return None

        return coordinates
    
    def calc_center_of_mass(coordinates):
        '''Calculates the center of mass of the molecule. Returns CoM along the x, y, and z axes'''

        #calculate the total mass of the molecule
        try:
            total_mass = sum(atom_masses[atom] for (atom, _, _, _) in coordinates)
        
        #if anything non-numeric gets passed, then sum will throw a value error
        except ValueError:
            return None
        
        #calculate the center of mass along the x, y, and z axes
        try:
            center_x = sum(atom_masses[atom] * float(x) for (atom, x, _, _) in coordinates) / total_mass
            center_y = sum(atom_masses[atom] * float(y) for (atom, _, y, _) in coordinates) / total_mass
            center_z = sum(atom_masses[atom] * float(z) for (atom, _, _, z) in coordinates) / total_mass
        
        #if anything non-numeric gets passed, then float/sum will throw a value error
        except ValueError:
            return None
        
        #have the function return the x, y, and z centers of mass
        return center_x, center_y, center_z

    def mass_weighted_distances(coordinates, center_of_mass):
        '''Calculates the distance of each atom in a molecule from the molecular center of mass, and weights that distance by atomic mass'''

        #initialize vector filled with zeroes with the same length as the number of atoms
        mw_distances = np.zeros(len(coordinates))
        
        #loop through each atom in the coordinates, calculate its distnace from the center of mass, then weight that distance according to the atomic weight
        for i, (atom, x, y, z) in enumerate(coordinates):
            distance = np.sqrt(np.square(float(x) - center_of_mass[0]) + np.square(float(y) - center_of_mass[1]) + np.square(float(z) - center_of_mass[2]))
            mw_distances[i] = atom_masses[atom] * np.abs(distance)

        #sort distnaces from smallest to largest. If we don't do this, files with identical geometries but different arrangements of atoms can be mistakenly assigned as unique.
        sorted_indices = np.argsort(mw_distances)
        sorted_distances = mw_distances[sorted_indices]
        
        #return the 1D vector of mass-weighted distances of each atom from the molecular center of mass
        return sorted_distances

    def compute_vector_similarity(vector1, vector2):
        '''Computes the similiarity between two vectors by the angle between them'''
        #compute dot product of two vectors, as well as their magnitudes
        dot_product = np.dot(vector1, vector2)
        magnitude1 = np.linalg.norm(vector1)
        magnitude2 = np.linalg.norm(vector2)
        
        #compute angle between the two vectors (this is the cosine similarity!)
        similarity = dot_product / (magnitude1 * magnitude2)

        #for details on this equation, please see the following publication: https://www.frontiersin.org/articles/10.3389/fchem.2019.00519/full
        cosine_similarity = (1. - ((1/np.pi) * np.arccos(similarity))) * 100

        #return the similiarity between the two vectors
        return cosine_similarity

    def main(directory_or_csv, similarity_threshold, write_to_file):
        '''Main workflow'''

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting cosine similiarity sorting of the files found in {os.path.basename(directory_or_csv)} ...')
        QApplication.processEvents()

        #Check if input is a CSV file or a directory
        if directory_or_csv.endswith('.csv'):
            destination = os.path.dirname(directory_or_csv)
            with open(directory_or_csv, 'r') as csvfile:
                reader = csv.reader(csvfile)
                next(reader)
                
                #get the file list from the .csv
                file_list = [row[0] for row in reader]

        else:
            destination = directory_or_csv
            
            #sort file list in natural order, not funky python string comprehension way
            file_list = natsorted([file for file in os.listdir(directory_or_csv) if file.lower().endswith('.gjf') or file.lower().endswith('.inp') or file.lower().endswith('.xyz')]) 

        #empty list to write distance vectors to
        vectors = []

        #all molecules should have the same number of atoms, but we'll check for this anyways
        expected_num_atoms = None

        #Start conversion of atomic coordinates to vectors
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Converting molecular coordinates to a 1D vector composed of mass-weighted atomic distances to the molecular center of mass, sorted from smallest to largest values...')
        QApplication.processEvents()   

        for file_name in file_list:

            file_path = os.path.join(destination, file_name)
            coordinates = read_file(file_path)
            
            #if coordinates were not extracted from the file, let the user know, and exit the processing
            if coordinates is None:
                
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No coordinates were found in {file_name}. Processing will now be aborted - please check the problematic file.')
                QApplication.processEvents()  
                return
            
            #Check the number of atoms in the current molecule, ensuring that it matches the previous file
            current_num_atoms = len(coordinates)
            
            if expected_num_atoms is None:
                #Initialize expected_num_atoms with the number of atoms in the first molecule
                expected_num_atoms = current_num_atoms
            
            elif current_num_atoms != expected_num_atoms:
                #Raise an error if the number of atoms is not consistent
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Number of atoms in {file_name} is inconsistent with the number of atoms in previous files. Please check that your .gjf files contain reasonable structures.')
                QApplication.processEvents()
                return

            CoM = calc_center_of_mass(coordinates)

            #throw an error if the CoM could not be calculated
            if CoM is None:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} A value error was encountered while calcualting the CoM of {file_name}. Processing will now be aborted - please check the problematic file for non-numeric characters in the coordinate block.')
                QApplication.processEvents()  
                return
                            
            vector = mass_weighted_distances(coordinates, CoM)
            vectors.append([file_name,vector])

        #determine the similarity between all vectors in the list by computing pairwise cosine similiarities. 

        #List to store indices of duplicate vectors and calculated pairwise similiarities
        duplicates = []
        similarity_list = []

        #Make a copy of vectors to preserve the original list
        unique_vectors = vectors

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Sorting vectors based on cosine similiarity...')
        QApplication.processEvents()   

        #Iterate through the sorted list of vectors
        for idx_i, (file_name_i, vector_i) in enumerate(vectors):
            #Skip if the vector is already marked as a duplicate
            if idx_i in duplicates:
                continue

            #Compare the current vector with others
            for idx_j, (file_name_j, vector_j) in enumerate(vectors):
                
                #Skip comparison with itself or if the vector is already marked as a duplicate
                if idx_j <= idx_i or idx_j in duplicates:
                    continue

                try:
                    #Calculate cosine similarity between vectors
                    similarity = compute_vector_similarity(vector_i, vector_j)
                except ValueError:
                    print(f'{file_name_i} and {file_name_j} are vectors of different sizes and cannot be compared. Are they isomers of one antoher? Check the file')
                    continue  #If structures of two different sizes are compared

                #Mark the vector as duplicate if similarity is above the threshold, and specify in the similiarity list that filename_j is not unique compared to filename_i
                if similarity > similarity_threshold:
                    duplicates.append(idx_j)
                    similarity_list.append([file_name_i, file_name_j, similarity, 'No'])
                
                else: 
                    similarity_list.append([file_name_i, file_name_j, similarity, 'Yes'])
                    
        #Remove duplicate vectors from the unique vectors list
        for index in sorted(duplicates, reverse=True):
            del unique_vectors[index]

        #Check if output directory exists, create if not
        base_output_directory = os.path.join(destination, f'uniques_sim_{str(similarity_threshold).replace(".", "-")}')
        output_directory = base_output_directory

        i = 2
        while os.path.exists(output_directory):
            output_directory = f'{base_output_directory}_v{i}'
            i += 1

        os.makedirs(output_directory)

        #copy unique files to unqiues directory
        if directory_or_csv.lower().endswith('.csv'):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(unique_vectors)} unique geometries identified by cosine similiarty and energetic sorting are being copied to {os.path.basename(output_directory)}:')

        else:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(unique_vectors)} unique geometries identified by cosine similiarty are being copied to {os.path.basename(output_directory)}:')

        for file, _ in unique_vectors:
            #print(file)
            shutil.copy(os.path.join(destination, file), os.path.join(output_directory, file))

        if write_to_file:
            sim_file = os.path.join(output_directory, 'pairwise_similarities.csv')

            #header for the .csv
            header = ['Molecule 1', 'Molecule 2', 'Cosine Similarity', 'Is molecule 2 unique?']
            
            #convert similiarity list to a numpy array and round the similiarities to 2 decimal places
            similarity_list = np.array(similarity_list)
            similarity_list[:, 2] = np.round(similarity_list[:, 2].astype(float), 5)
            
            #Concatenate the header as the first row
            similarity_list_with_header = np.vstack([header, similarity_list])
            
            #Save the array to a text file with consistent spacing
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Writing pairwise similiarities to {os.path.join(os.path.basename(output_directory), os.path.basename(sim_file))}.')
            np.savetxt(sim_file, similarity_list_with_header, fmt='%-30s,%-30s,%-30s,%-30s')

    #execute main function to calculate molecular similarity
    stime = time.time()
    main(directory_or_csv, similarity_threshold, write_to_file)    
    ftime = time.time()

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Similarity analysis was completed in {np.round((ftime - stime),2)} seconds.')
    return

#for external testing
if __name__ == "__main__":
    directory_or_csv = r'C:\Users\Scott Hopkins\Downloads\Sample_Files\T2_Cosine_sim\Conformers_gjf\Energies.csv'
    similarity_threshold = 98.5
    write_to_file = True

    #to run as a standalone, you need to change the location of the atom masses file
    #change: from Python.atom_mass import atom_masses
    #to: from atom_mass import atom_masses
    cosine_sim(directory_or_csv, similarity_threshold, write_to_file)
    