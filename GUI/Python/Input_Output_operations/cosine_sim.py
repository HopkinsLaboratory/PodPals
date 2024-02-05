#import os, time, csv, shutil
import numpy as np
from natsort import natsorted
from datetime import datetime
from Python.atom_mass import atom_masses

def cosine_sim(directory_or_csv, similarity_threshold, write_to_file):
    
    '''
    Calculate cosine similarity between molecular structures in a directory or specified in a CSV file.

    Parameters:
    - directory_or_csv: Path to a directory containing .gjf files or a CSV file with a list of file paths.
    - similarity_threshold: Threshold for considering molecules as similar.
    - write_to_file: Boolean flag to indicate whether to write results to a file.
    '''

    def read_gjf(file_path):
        """Read the xyz data of a .gjf file and write the data to a list."""
        
        #empty list to write coordinate data to
        coordinates = []

        with open(file_path, 'r') as file:
            content = file.readlines()

        # Find the index of the line containing two integers separated by a space
        cutoff_index = next((i for i, line in enumerate(content) if len(line.split()) == 2), None)

        # If cutoff line is found, return content starting from the next line
        if cutoff_index is not None:
            xyz_data = content[cutoff_index + 1:]
        
            for line in xyz_data:
                line = line.split()
                
                if len(line) == 4: #only lines containing the atom type and x, y, and z coordinates will have a length of 4
                    coordinates.append(line)
                
            return coordinates
        
        else:
            # Handle the case where no cutoff line is found
            raise ValueError(f'The charge and multiplicity was not found in {os.path.basename(file_path)}. Please check that this file is in the .gjf format.')

    def calc_center_of_mass(coordinates):
        '''Calculates the center of mass of the molecule. Returns CoM along the x, y, and z axes'''

        #calculate the total mass of the molecule
        total_mass = sum(atom_masses[atom] for (atom, _, _, _) in coordinates)
        
        #calculate the center of mass along the x, y, and z axes
        center_x = sum(atom_masses[atom] * float(x) for (atom, x, _, _) in coordinates) / total_mass
        center_y = sum(atom_masses[atom] * float(y) for (atom, _, y, _) in coordinates) / total_mass
        center_z = sum(atom_masses[atom] * float(z) for (atom, _, _, z) in coordinates) / total_mass
        
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

        # sort distnaces from smallest to largest. If we don't do this, files with identical geometries but different arrangements of atoms can be mistakenly assigned as unique.
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

        #scale the similiarity such that it returns a number understandable by our human brains. 
        #With this, two vectors that are identical return a value of 100 (100%) similiar. Perpindicular vectors yield a value of 0 (0% similiar), and equal but opposite vectors yield a value of -100 (totally opposite).
        similarity_scale = (200/np.pi) * np.arcsin(similarity)
        similarity_scale = (1. - ((1/np.pi) * np.arccos(similarity))) * 100

        #return the similiarity between the two vectors
        return similarity_scale

    def main(directory_or_csv, similarity_threshold, write_to_file):
        '''Main workflow'''
        
        # Check if input is a CSV file or a directory
        if directory_or_csv.endswith('.csv'):
            destination = os.path.dirname(directory_or_csv)
            with open(directory_or_csv, 'r') as csvfile:
                reader = csv.reader(csvfile)
                
                #read the first line in the .csv
                first_line = next(reader)
                
                #look at the first line in the .csv and check if its a header or a list of filenames
                #entry1 = first_line.split()[0]
                
                if first_line[0].lower().endswith('.gjf'):
                    file_list = [first_line[0]] + [row[0] for row in reader]

                else:
                    file_list = [row[0] for row in reader]

        else:
            destination = directory_or_csv
            file_list = natsorted([file for file in os.listdir(directory_or_csv) if file.lower().endswith('.gjf')]) #sort file list in natural order, not funky python string comprehension way

        #empty list to write distance vectors to
        vectors = []

        # all molecules should have the same number of atoms, but we'll check for this anyways
        expected_num_atoms = None

        for file_name in file_list:

            file_path = os.path.join(destination, file_name)
            coordinates = read_gjf(file_path)
            # Check the number of atoms in the current molecule
            current_num_atoms = len(coordinates)
            
            if expected_num_atoms is None:
                # Initialize expected_num_atoms with the number of atoms in the first molecule
                expected_num_atoms = current_num_atoms
            elif current_num_atoms != expected_num_atoms:
                # Raise an error if the number of atoms is not consistent
                raise ValueError(f'Number of atoms in {file_name} is inconsistent with the number of atoms in previous .gjf files.')

            CoM = calc_center_of_mass(coordinates)
            vector = mass_weighted_distances(coordinates, CoM)
            vectors.append([file_name,vector])

        #determine the similarity between all vectors in the list by computing pairwise cosine similiarities. 

        # List to store indices of duplicate vectors and calculated pairwise similiarities
        duplicates = []
        similarity_list = []

        # Make a copy of vectors to preserve the original list
        unique_vectors = vectors

        # Iterate through the sorted list of vectors
        for idx_i, (file_name_i, vector_i) in enumerate(vectors):
            # Skip if the vector is already marked as a duplicate
            if idx_i in duplicates:
                continue

            # Compare the current vector with others
            for idx_j, (file_name_j, vector_j) in enumerate(vectors):
                
                # Skip comparison with itself or if the vector is already marked as a duplicate
                if idx_j <= idx_i or idx_j in duplicates:
                    continue

                try:
                    # Calculate cosine similarity between vectors
                    similarity = compute_vector_similarity(vector_i, vector_j)
                    similarity_list.append([file_name_i, file_name_j, similarity])
                except ValueError:
                    print(f'{file_name_i} and {file_name_j} are vectors of different sizes and cannot be compared. Are they isomers of one antoher? Check the .gjf file')
                    continue  # If structures of two different sizes are compared

                # Mark the vector as duplicate if similarity is above the threshold
                if similarity > similarity_threshold:
                    duplicates.append(idx_j)

        # Remove duplicate vectors from the unique vectors list
        for index in sorted(duplicates, reverse=True):
            del unique_vectors[index]

        # Check if output directory exists, create if not
        output_directory = os.path.join(destination, 'uniques')
        i = 1
        while os.path.exists(output_directory):
            try:
                output_directory = f'{os.path.join(destination, f"uniques_{i}")}' # Alex
                os.makedirs(output_directory)
                break  # Break the loop if the directory is successfully created
            except FileExistsError:
                i += 1

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        # copy unique files to unqiues directory
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(unique_vectors)} unique geometries based on cosine similarity and energy have been copied to {os.path.basename(output_directory)}:')
        for file, _ in unique_vectors:
            #print(file)
            shutil.copy(os.path.join(destination, file), os.path.join(output_directory, file))

        if write_to_file:
            sim_file = os.path.join(output_directory, 'similarities.csv')

            #header for the .csv
            header = ['Molecule 1', 'Molecule 2', 'Cosine Similarity']
            
            #convert similiarity list to a numpy array and round the similiarities to 2 decimal places
            similarity_list = np.array(similarity_list)
            similarity_list[:, 2] = np.round(similarity_list[:, 2].astype(float), 5)
            
            # Concatenate the header as the first row
            similarity_list_with_header = np.vstack([header, similarity_list])
            
            # Save the array to a text file with consistent spacing
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Writing pairwise similiarities to {os.path.join(os.path.basename(output_directory), os.path.basename(sim_file))}.')
            np.savetxt(sim_file, similarity_list_with_header, fmt='%-30s,%-30s,%-30s')

    #execute main function to calculate molecular similarity
    stime = time.time()
    main(directory_or_csv, similarity_threshold, write_to_file)    
    ftime = time.time()

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Similarity analysis was completed in {np.round((ftime - stime),2)} seconds.')

# for external testing
if __name__ == "__main__":
    directory_or_csv = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\MobCal-MPI\MobCal-MPI\Manual\Appendix_A\CREST_Outputs\cosine_sim_testing_3\Energies.csv'
    similarity_threshold = 99.
    write_to_file = True

    #to run as a standalone, you need to change the location of the atom masses file
    #change: from Python.atom_mass import atom_masses
    #to: from atom_mass import atom_masses
    cosine_sim(directory_or_csv, similarity_threshold, write_to_file)
    