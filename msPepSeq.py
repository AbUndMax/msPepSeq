# Program for Peptide sequencing based on b- and y-ion fragmentation
# Author: Niklas G.

# USAGE:
# python msPepSeq.py <path-to-ms-spectra> [--output_to_file] [--output_assignments]

# the last two flags are optionals that enable:
# to output the candidate sequences to a file in .out/ instead to the console
# to output all assignemnts that were calculated per peak for b and y ions to respective files in .out/




from scipy.spatial.distance import euclidean
import argparse
import os

# Monoisotopic masses for amino acids (One-Letter Code) from the lecture
AA_TO_MASS_LECTURE = {
    'A': 71.03711,   
    'R': 156.10111,  
    'N': 114.04293,  
    'D': 115.02694,  
    'C': 103.00919,  
    'E': 129.04259,  
    'Q': 128.05858,  
    'G': 57.02146,   
    'H': 137.05891,  
    'I': 113.08406,  
    'L': 113.08406,  
    'K': 128.09496,  
    'M': 131.04049,  
    'F': 147.06841,  
    'P': 97.05276,   
    'S': 87.03203,   
    'T': 101.04768,  
    'W': 186.07931,  
    'Y': 163.06333,  
    'V': 99.06841    
}

# Define mass tolerance as written in assignment
MASS_TOLERANCE = 0.055  # Da
WATER = 18.01056

class Peak:
    def __init__(self, m_z: float):
        self.m_z = m_z
        self.b_mass: float = m_z - 1  # (M + 1 in MS)
        self.y_mass: float = m_z - 19 # (M + 19 in MS because from H added to NH2, OH on CO and a proton when fragmented)
        self.possible_b_peptides: list[str] = []
        self.possible_y_peptides: list[str] = []
        
        
def load_peaks(filePath):
    """
    Load m/z values from a given file, create Peak objects for each, 
    and return a list of peaks including a synthetic start peak.

    Parameters:
    filePath (str): Path to the input file containing m/z values.

    Returns:
    tuple: A list of Peak objects (including start peak) and a list of raw m/z values.
    """
    # Read m/z values from the input file and create Peak objects
    peaks = []
    masses = []
    with open("../materials_08/b_y_spectrum.txt") as f:
        for line in f:
            # Each line contains a single m/z value; strip whitespace and convert to float
            m_z = float(line.strip())
            masses.append(m_z)
            peaks.append(Peak(m_z))
            
    # Create synthetic start peak (b_mass = 0) to anchor sequence building
    start = Peak(m_z=0.0)
    start.b_mass = 0
    start.y_mass = 0
    start.possible_b_peptides = [("", "source_peak")] # the second entry in tuple is just to see where the peptide came from
    start.possible_y_peptides = [("", "source_peak")]

    # Prepend the start peak and sort all peaks by ascending b-ion mass
    # Sorting ascendingly ensures we always compare each peak only against lower-mass ones.
    return [start] + sorted(peaks, key=lambda peak: peak.b_mass), sorted(masses)


def assign_peptides_to_peaks(peaks):
    """
    Assign possible peptide sequences to each peak based on pairwise mass differences.
    Matches amino acid masses to differences within specified tolerance.

    Parameters:
    peaks (list of Peak): List of Peak objects to process.
    """
    # using "pairwise difference analysis" -> take peak, compare it with all other peaks below
    for i in range(1, len(peaks)):
        peak = peaks[i]
        
        for j in range(0, i):
            peak_below = peaks[j]
            
            ########################################################## b-ion processing
            # Calculate mass difference for b-ions
            diff_b = peak.b_mass - peak_below.b_mass
            
            # Find all AAs whose mass matches this diff within tolerance
            aas_to_add_to_b = []
            for aa, aa_mass in AA_TO_MASS_LECTURE.items():
                if abs(diff_b - aa_mass) <= MASS_TOLERANCE * 2: 
                    # *2 because we compare two peaks. tolerance is applied on each peak
                    # The tolerance times two makes sense since in the worst case the peak with bigger mass could have at the upper end of
                    # the tolerance while the lower peak mass could have been at the lower end of the tolereance and thus by substraction
                    # will result in a bigger difference than any aa with a single tolerance could achieve
                    aas_to_add_to_b.append(aa)
            
            # Extend each candidate b-peptide from the lower peak by those AAs
            if aas_to_add_to_b:
                for aa in aas_to_add_to_b:
                    for peptide, _ in peak_below.possible_b_peptides:
                        peak.possible_b_peptides.append((peptide + aa, f"from: {peak_below.b_mass} diff: {diff_b}"))
            
            ########################################################## y-ion processing  
            # Calculate mass difference for y-ions
            diff_y = peak.y_mass - peak_below.y_mass
            
            # Find all AAs whose mass matches this diff within tolerance
            aas_to_add_to_y = []
            for aa, aa_mass in AA_TO_MASS_LECTURE.items():
                if abs(diff_y - aa_mass) <= MASS_TOLERANCE * 2:
                    aas_to_add_to_y.append(aa)

            # Extend each candidate b-peptide from the lower peak by those AAs
            if aas_to_add_to_y:
                for aa in aas_to_add_to_y:
                    for peptide, _ in peak_below.possible_y_peptides:
                        peak.possible_y_peptides.append((aa + peptide, f"from: {peak_below.y_mass} diff: {diff_y}"))
                        
                        
def save_assignments_to_file(peaks):
    """
    Save the assigned possible peptides for each peak to text files for b-ions and y-ions.

    Parameters:
    peaks (list of Peak): List of Peak objects containing peptide assignments.
    """
    with open(os.path.join('out', "out_b.txt"), "w") as f:
        f.write("possible peptides:\n\n")
        for peak in peaks:
            f.write(f"peak: {peak.m_z} in b: {peak.b_mass}\n")
            for peptide in peak.possible_b_peptides:
                f.write(f"pep: {peptide}\n")
            f.write("\n")
            
    with open(os.path.join('out', "out_y.txt"), "w") as f:
        f.write("possible peptides:\n\n")
        for peak in peaks:
            f.write(f"peak: {peak.m_z} in y: {peak.y_mass}\n")
            for peptide in peak.possible_y_peptides:
                f.write(f"pep: {peptide}\n")
            f.write("\n")


def find_common_sequences(peaks, peptide_length):
    """
    Find peptide sequences that are common to both b-ion and y-ion paths 
    with the expected peptide length.

    Parameters:
    peaks (list of Peak): List of Peak objects with assigned peptide sequences.
    peptide_length (int): Expected length of the peptide sequence.

    Returns:
    list: Candidate peptide sequences that match both b- and y-ion assignments.
    """
    # find common sequences:
    candidates = []

    # the last peak is the longest y-ion
    for y_peptide, _ in peaks[-1].possible_y_peptides:
        
        # the peak before the last peak is the longest b-ion (it is exactly 18Da away from the longest y-ion (differnece of water))
        for b_peptide, _ in peaks[-2].possible_b_peptides:
            
            if y_peptide == b_peptide and len(y_peptide) == peptide_length:
                candidates.append(y_peptide)
    
    return candidates


def report_candidates(filePath, candidates, peptide_length):
    """
    Print or save the list of candidate peptide sequences with their distances.

    Parameters:
    filePath (str): If provided, save the output to a file; otherwise, print to console.
    candidates (list of tuple): Candidate peptide sequences and their distances.
    peptide_length (int): Expected length of the peptide sequence.
    """
    output = f"> Found {len(candidates)} candidates of length: {peptide_length}\n"
    output += "> For similarity measure, euclidean distance based scoring was used (1 / (1 + eucl_dist))\n"
    
    for candidate, distance in candidates:
        output += candidate + " - score: " + str(distance) + "\n"
    
    if filePath:
        with open(os.path.join('out', "candidate_report.txt"), 'w') as f:
            f.write(output)
    else:
        print(output)
        
        
        
        
        
# Similarity assessment:
def calculate_theoretical_ms(peptide):
    """
    Calculate the theoretical m/z values for b- and y-ions of a given peptide.

    Parameters:
    peptide (str): Peptide sequence.

    Returns:
    list: Sorted list of theoretical b- and y-ion masses.
    """
    def calculate_mass(peptide):
        mass = 0
        for aa in peptide:
            mass += AA_TO_MASS_LECTURE[aa]    
        return mass
    
    masses = []
    for i in range(len(peptide) + 1):
        b_fragment = peptide[:i]
        y_fragment = peptide[i:len(peptide)]
        
        b_mass = calculate_mass(b_fragment) + 1
        y_mass = calculate_mass(y_fragment) + 19
        
        # exlcude the "starter" peptides but append all others
        if b_mass != 1: masses.append(b_mass)
        if y_mass != 19: masses.append(y_mass)
        
    return sorted(masses)


def assess_distances_to_input(candidates, input_ms):
    """
    Calculate the Euclidean distance between theoretical and observed m/z spectra for each candidate.

    Parameters:
    candidates (list of str): List of candidate peptide sequences.
    input_ms (list of float): List of observed m/z values.

    Returns:
    list: List of tuples containing peptide sequences and their distances, sorted descending by distance.
    """
    distances = []
    for candidate in candidates:
        dist = euclidean(calculate_theoretical_ms(candidate), input_ms)
        score = 1 / (1 + dist)
        distances.append((candidate, score))
        
    distances.sort(key=lambda x: x[1], reverse=True)
    return distances



def main():
    parser = argparse.ArgumentParser(description="Process b/y spectrum peak positions from a text file.")
    parser.add_argument("input_file", type=str, help="Path to a txt file containing one m/z value per line")
    parser.add_argument("--output_to_file", action="store_true", help="if specified, output is written to file")
    parser.add_argument("--output_assignments", action="store_true", help="saves the calculated assigned possible peptides per peak to files")
    args = parser.parse_args()
    
    # ensure out folder exists if needed
    if args.output_to_file or args.output_assignments:
        os.makedirs("out", exist_ok=True)
    
    # load peaks
    peaks, input_ms = load_peaks(args.input_file)
    peptide_length = len(peaks) // 2 # just divide by two because we have exactly the same amount of b and y ions
    
    # assigne peptides per peak
    assign_peptides_to_peaks(peaks)
    
    # save assignemnt if specified:
    if args.output_assignments: save_assignments_to_file(peaks)
    
    # calculate candidates
    candidates = find_common_sequences(peaks, peptide_length)
    
    # similarity assessment - how similar is the ms of the constructed sequence to the input ms:
    candidates = assess_distances_to_input(candidates, input_ms)
    
    # report candidates
    report_candidates(args.output_to_file, candidates, peptide_length)
    

if __name__ == "__main__":
    main()