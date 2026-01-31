#!/usr/bin/env python3
"""
Extract sphere centers from CSV file and generate OBJ file for specific spheres.
"""

import csv
import sys

def extract_sphere_centers(csv_file, sphere_indices, output_obj):
    """
    Extract sphere centers for specified indices and create OBJ file.
    
    Args:
        csv_file: Path to input CSV file containing sphere data
        sphere_indices: List of sphere indices to extract
        output_obj: Path to output OBJ file
    """
    sphere_centers = {}
    
    # Read CSV file and extract sphere data
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        
        # Skip header if present
        first_row = next(reader, None)
        if first_row and first_row[0].lower() in ['x', 'x1coord', 'xcoord']:
            # Header row, continue with data
            pass
        else:
            # No header, rewind and process first row as data
            f.seek(0)
            reader = csv.reader(f)
        
        # Read sphere data
        for row_num, row in enumerate(reader):
            if len(row) >= 4:  # x, y, z, r
                sphere_idx = row_num
                if sphere_idx in sphere_indices:
                    x, y, z, r = float(row[0]), float(row[1]), float(row[2]), float(row[3])
                    sphere_centers[sphere_idx] = (x, y, z, r)
    
    # Write OBJ file
    with open(output_obj, 'w') as f:
        f.write("# OBJ file with sphere centers\n")
        f.write(f"# Spheres: {sphere_indices}\n")
        f.write("# Format: v x y z\n")
        
        for idx in sorted(sphere_indices):
            if idx in sphere_centers:
                x, y, z, r = sphere_centers[idx]
                f.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
                print(f"Sphere {idx}: center=({x:.6f}, {y:.6f}, {z:.6f}), radius={r:.6f}")
            else:
                print(f"Warning: Sphere {idx} not found in CSV file")
    
    print(f"\nGenerated OBJ file: {output_obj}")
    print(f"Total spheres extracted: {len(sphere_centers)}")

def main():
    # Configuration
    csv_file = "/home/yiming/research/vorocrust/data/spheres/Sphere_2500_55.csv"
    sphere_indices = [0, 2, 3, 723]
    output_obj = "sphere_centers_0_2_3_723.obj"
    
    # Extract sphere centers and generate OBJ
    extract_sphere_centers(csv_file, sphere_indices, output_obj)

if __name__ == "__main__":
    main()
