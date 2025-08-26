#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced Fixed Trajectory Analysis to Interactive HTML Generator
Main entry point - coordinates all modules
"""

import argparse
import sys
from pathlib import Path

# Import from parallel modules
from contact_analyzer import FastContactAnalyzer
from html_builder import HTMLBuilder
from data_processor import DataProcessor
from visualization_generator import VisualizationGenerator
import mdtraj as md

def main():
    parser = argparse.ArgumentParser(description="Enhanced trajectory analysis to interactive HTML")
    parser.add_argument("trajectory", help="Trajectory file (.xtc, .dcd, etc.)")
    parser.add_argument("topology", help="Topology file (.pdb, .gro, etc.)")
    parser.add_argument("ligand", help="Ligand file (.sdf, .mol, .mol2)")
    parser.add_argument("-o", "--output", default="contact_analysis.html", help="Output HTML file")
    
    args = parser.parse_args()
    
    print("=== Enhanced Contact Analysis to HTML ===")
    
    try:
        # Load trajectory
        print(f"Loading trajectory: {Path(args.trajectory).name}")
        traj = md.load(args.trajectory, top=args.topology)
        print(f"Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")
        
        # Load and process ligand
        data_processor = DataProcessor()
        ligand_mol = data_processor.load_ligand_structure(args.ligand)
        
        # Analyze contacts
        print("\nAnalyzing contacts...")
        analyzer = FastContactAnalyzer(traj)
        contact_results = analyzer.calculate_contact_proportions()
        
        # Generate visualization data
        vis_generator = VisualizationGenerator()
        ligand_data = vis_generator.generate_ligand_data(contact_results, ligand_mol)
        contacts = vis_generator.generate_contact_data(contact_results, ligand_data)
        
        # Calculate statistics
        stats = data_processor.calculate_statistics(contacts)
        
        # Generate HTML
        print("\nGenerating interactive HTML...")
        html_builder = HTMLBuilder()
        html_content = html_builder.generate_html(
            trajectory_file=Path(args.trajectory).name,
            topology_file=Path(args.topology).name,
            ligand_file=Path(args.ligand).name,
            ligand_name=contact_results['ligand_residue'].name if contact_results['ligand_residue'] else 'UNK',
            total_frames=contact_results['total_frames'],
            ligand_data=ligand_data,
            contacts=contacts,
            stats=stats
        )
        
        # Write output
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"Interactive HTML generated: {args.output}")
        print(f"Found {stats['total_contacts']} significant contacts")
        print(f"Open the file in a web browser to view and adjust the layout")
        print(f"\n=== Complete! ===")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()