import os
import sys

# Parse command-line arguments
diamond_file = sys.argv[1]   # DIAMOND alignment output
faa_file = sys.argv[2]       # Protein sequence file in FASTA format
result_file = sys.argv[3]    # Output result file

# Dictionary to store protein counts per genome
genome_protein_counts = {}

# First pass: Count proteins per genome from FASTA headers
with open(faa_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            # Extract genome name by removing the last suffix (e.g. >genome_1 -> genome)
            protein_id = line[1:].split()[0]  # Get first field after '>'
            genome_name = '_'.join(protein_id.split('_')[:-1])
            
            # Update protein count for this genome
            genome_protein_counts[genome_name] = genome_protein_counts.get(genome_name, 0) + 1

# Dictionary to store best alignment scores
# Structure: {query_genome: {target_genome: {query_protein: best_score}}}
alignment_data = {}

# Second pass: Initialize alignment_data structure with genomes
with open(faa_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            protein_id = line[1:].split()[0]
            genome_name = '_'.join(protein_id.split('_')[:-1])
            if genome_name not in alignment_data:
                alignment_data[genome_name] = {}

# Process DIAMOND alignment results
processed_count = 0
with open(diamond_file, 'r') as f:
    for line in f:
        processed_count += 1
        if processed_count % 100000 == 0:
            print(f"Processed {processed_count} alignments")
        
        fields = line.strip().split('\t')
        query_protein = fields[0]
        target_protein = fields[1]
        identity = float(fields[2])  # Alignment identity percentage
        
        # Extract genome names from protein IDs
        query_genome = '_'.join(query_protein.split('_')[:-1])
        target_genome = '_'.join(target_protein.split('_')[:-1])
        
        # Skip if query genome not in our database
        if query_genome not in alignment_data:
            continue
            
        # Initialize nested dictionaries if missing
        if target_genome not in alignment_data[query_genome]:
            alignment_data[query_genome][target_genome] = {}
        
        # Convert identity to fraction (0.0-1.0 range)
        normalized_score = identity / 100.0
        
        # Track best score per query protein
        current_best = alignment_data[query_genome][target_genome].get(query_protein, -1)
        if normalized_score > current_best:
            alignment_data[query_genome][target_genome][query_protein] = normalized_score

# Write final results
with open(result_file, 'w') as out:
    # Header for output file
    out.write("Query_Genome\tTarget_Genome\tQuery_Proteins\tTarget_Proteins\tAligned_Proteins\tNormalized_Score\n")
    
    for query_genome, targets in alignment_data.items():
        for target_genome, proteins in targets.items():
            # Skip self-comparisons and invalid pairs
            if query_genome == target_genome:
                continue
            if query_genome not in genome_protein_counts or target_genome not in genome_protein_counts:
                continue
            
            # Calculate metrics
            total_score = sum(proteins.values())
            aligned_count = len(proteins)
            min_proteins = min(genome_protein_counts[query_genome], genome_protein_counts[target_genome])
            normalized_score = total_score / min_proteins
            
            # Write results
            out.write(f"{query_genome}\t{target_genome}\t")
            out.write(f"{genome_protein_counts[query_genome]}\t{genome_protein_counts[target_genome]}\t")
            out.write(f"{aligned_count}\t{normalized_score:.6f}\n")
